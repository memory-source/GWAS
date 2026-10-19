"""Governance prototype: reveal manifest + egress monitor + hash-chained ledger.

This makes the release-object layer of the governance mechanism concrete (the
paper otherwise argues it, not builds it). It is a PLAINTEXT, release-object
level prototype:

  - a declarative, content-hashed, versioned manifest of what a run is allowed
    to release (per stripped quantity: semantic type, precision, recipient,
    condition);
  - an EgressMonitor that wraps a `release_objects()` dict and can only emit
    quantities declared in the manifest, coerced to the declared precision
    (decimal rounding, top-k truncation, binning, threshold-bool);
  - a tamper-evident (hash-chained) append-only ledger recording a VALUE HASH
    per release (never raw values), with the four cross-run relations from the
    proposal's Definition 3 (refinement / condition-loosening / duplicate /
    recipient-expansion).

It is NOT an MPC implementation and NOT "cryptographically signed" (no PKI /
recipient key custody): the manifest is content-hashed and versioned, and the
ledger is hash-chained. The sole-egress binding over this monitor is built in
`code/enforce.py` (a process-isolation runtime with a PEP-578 audit hook and a
single sanctioned pipe channel), tested against the careless-author bypass.
"""
from __future__ import annotations

import hashlib
import json
import os
import time
import numpy as np


def _canonical(obj):
    return json.dumps(obj, sort_keys=True, separators=(",", ":"))


def _sha(obj_or_bytes):
    b = obj_or_bytes if isinstance(obj_or_bytes, bytes) else _canonical(obj_or_bytes).encode()
    return "sha256:" + hashlib.sha256(b).hexdigest()


def _value_hash(arr):
    """Hash of a released array's canonical bytes (no raw values stored). A
    per-group dict is hashed entrywise in sorted-key order."""
    if isinstance(arr, dict):
        return _sha({str(k): _value_hash(v) for k, v in sorted(arr.items())})
    a = np.asarray(arr, dtype=np.float64)
    return _sha(a.tobytes() + _canonical(list(a.shape)).encode())


DEFAULT_QUANTITIES = ("af", "af_by_group", "lambda_gc", "eig", "loadings",
                      "missing", "pvals", "condition_bit")


def make_manifest(policy, run_id="run0", version=1):
    """policy: {quantity_id: {released: bool, precision: spec, recipient: str,
    condition: str}}. Unlisted quantities default to not-released."""
    entries = {}
    for q in DEFAULT_QUANTITIES:
        spec = dict(policy.get(q, {}))
        entries[q] = dict(
            released=bool(spec.get("released", False)),
            semantic_type=spec.get("semantic_type", q),
            precision=spec.get("precision", "full"),
            recipient=spec.get("recipient", "public"),
            condition=spec.get("condition", "always"),
        )
    manifest = dict(run_id=run_id, version=version, quantities=entries,
                    created=time.strftime("%Y-%m-%dT%H:%M:%S"))
    manifest["hash"] = _sha({k: v for k, v in manifest.items() if k != "hash"})
    return manifest


def apply_precision(value, spec):
    """Coerce a released value to the declared precision spec. A per-group
    mapping (dict of group -> array) is coerced entrywise, so the monitor can
    emit the per-site allele-frequency object through the same path."""
    if isinstance(value, dict) and spec not in (None, "full", "withheld"):
        return {k: apply_precision(v, spec) for k, v in value.items()}
    if spec in (None, "full"):
        return value
    if spec == "withheld":
        return None
    if isinstance(spec, (int, float)):            # decimals
        return np.round(np.asarray(value), int(spec))
    if isinstance(spec, dict):
        kind = spec.get("kind")
        if kind == "topk":
            v = np.asarray(value)
            k = int(spec["k"])
            flat = v.ravel()
            # default keeps the largest-magnitude entries (eigenvalues); a
            # p-value table instead keeps the SMALLEST p-values.
            if spec.get("order") == "small":
                order = np.argsort(flat)[:k]
            else:
                order = np.argsort(-np.abs(flat))[:k]
            out = np.full_like(flat, np.nan)
            out[order] = flat[order]
            return out.reshape(v.shape)
        if kind == "decimals":
            return np.round(np.asarray(value), int(spec["d"]))
        if kind == "bin":
            edges = np.asarray(spec["edges"], float)
            return np.digitize(np.asarray(value), edges)
        if kind == "bool":
            return (np.asarray(value) >= float(spec["threshold"])).astype(int)
    raise ValueError(f"unknown precision spec {spec!r}")


class Violation(Exception):
    pass


class EgressMonitor:
    def __init__(self, manifest, ledger):
        self.manifest = manifest
        self.ledger = ledger
        self.log = []

    def release(self, rel, quantity_id, recipient=None, condition=None):
        spec = self.manifest["quantities"].get(quantity_id)
        if spec is None or not spec["released"]:
            self.log.append(dict(quantity_id=quantity_id, outcome="blocked",
                                 reason="undeclared"))
            self.ledger.append(quantity_id, "blocked", self.manifest["run_id"])
            raise Violation(f"undeclared reveal: {quantity_id}")
        if recipient is not None and recipient != spec["recipient"]:
            self.log.append(dict(quantity_id=quantity_id, outcome="blocked",
                                 reason="recipient-expansion"))
            self.ledger.append(quantity_id, "blocked", self.manifest["run_id"])
            raise Violation(f"recipient expansion for {quantity_id}")
        if condition is not None and condition != spec["condition"]:
            self.log.append(dict(quantity_id=quantity_id, outcome="blocked",
                                 reason="condition-loosened"))
            self.ledger.append(quantity_id, "blocked", self.manifest["run_id"])
            raise Violation(f"condition loosened for {quantity_id}")
        raw = rel.get(quantity_id)
        if raw is None:
            self.log.append(dict(quantity_id=quantity_id, outcome="blocked",
                                 reason="absent"))
            raise Violation(f"quantity absent: {quantity_id}")
        declared = spec["precision"]
        got = apply_precision(raw, declared)
        # over-precision detection: if the declared precision is coarser than the
        # value the producer tried to emit, record a coercion
        outcome = "coarsened" if declared != "full" else "released"
        self.log.append(dict(quantity_id=quantity_id, outcome=outcome,
                             declared=declared))
        self.ledger.append(quantity_id, outcome, self.manifest["run_id"],
                           value_hash=_value_hash(got))
        return got


class Ledger:
    """Hash-chained append-only ledger. Each record links to the previous hash;
    records store a value HASH, never raw released values."""

    def __init__(self, path=None):
        self.path = path
        self.records = []
        self._prev = "sha256:genesis"

    def append(self, quantity_id, outcome, run_id, value_hash=None, meta=None):
        rec = dict(idx=len(self.records), quantity_id=quantity_id, outcome=outcome,
                   run_id=run_id, value_hash=value_hash, prev=self._prev,
                   ts=time.time(), meta=meta or {})
        rec["hash"] = _sha({k: v for k, v in rec.items() if k != "hash"})
        self._prev = rec["hash"]
        self.records.append(rec)
        if self.path:
            with open(self.path, "a") as f:
                f.write(_canonical(rec) + "\n")
        return rec

    def verify_chain(self):
        prev = "sha256:genesis"
        for rec in self.records:
            body = {k: v for k, v in rec.items() if k != "hash"}
            if rec["prev"] != prev or _sha(body) != rec["hash"]:
                return False
            prev = rec["hash"]
        return True

    @staticmethod
    def relations(records):
        """The proposal's Definition-3 cross-run relations among released
        quantities: duplicate / refinement / recipient-expansion /
        condition-loosening. Here detected from record metadata."""
        flags = []
        seen = {}
        for r in records:
            if r["outcome"] not in ("released", "coarsened"):
                continue
            key = (r["run_id"], r["quantity_id"])
            if key in seen:
                flags.append(f"duplicate({r['quantity_id']}@{r['run_id']})")
            seen[key] = r
            m = r.get("meta") or {}
            if m.get("refines"):
                flags.append(f"refinement({r['quantity_id']}:{m['refines']})")
            if m.get("recipient_expansion"):
                flags.append(f"recipient_expansion({r['quantity_id']})")
            if m.get("condition_loosened"):
                flags.append(f"condition_loosening({r['quantity_id']})")
        return flags


POLICIES = {
    # ungoverned: everything at full precision (what the paper's frontier attacks)
    "ungoverned": {q: dict(released=True, precision="full") for q in DEFAULT_QUANTITIES},
    # governed: canonical safe release (AF @1dp public; spectrum eigenvalues only;
    # loadings/af_by_group/missingness/pvals withheld)
    "governed": {
        "af": dict(released=True, precision=1, recipient="public"),
        "eig": dict(released=True, precision={"kind": "topk", "k": 10}),
        "lambda_gc": dict(released=True, precision=2),
        "loadings": dict(released=False),
        "af_by_group": dict(released=False),
        "missing": dict(released=False),
        "pvals": dict(released=False),
    },
}


def violation_injection_suite():
    """Attempt the P0-5 A2 injections against the governed manifest and classify
    each outcome. Returns the matrix + a TCB inventory."""
    led = Ledger()
    mon = EgressMonitor(make_manifest(POLICIES["governed"], "inj"), led)
    # a fake release dict carrying everything a producer might try to emit
    rel = dict(af=np.array([0.51, 0.49]), eig=np.arange(10.0), lambda_gc=1.02,
               loadings=np.random.default_rng(0).normal(size=(2, 3)),
               af_by_group={0: np.array([0.5])}, missing=np.array([0.01]),
               pvals=np.array([1e-8]), condition_bit=1)
    matrix = []

    def attempt(name, qid, **kw):
        try:
            mon.release(rel, qid, **kw)
            matrix.append(dict(injection=name, outcome="released"))
        except Violation as e:
            matrix.append(dict(injection=name, outcome="blocked", reason=str(e)))

    attempt("undeclared_loadings", "loadings")
    attempt("undeclared_missingness", "missing")
    attempt("undeclared_condition_bit", "condition_bit")
    attempt("recipient_expansion", "af", recipient="adversary")
    attempt("condition_loosened", "af", condition="always+")
    # over-precision: producer emits full AF but manifest says 1dp -> coarsened
    got = mon.release(rel, "af")
    matrix.append(dict(injection="af_overprecision_to_full",
                       outcome="coarsened", emitted_precision=int(1),
                       returned_matches_declared=bool(
                           np.allclose(got, np.round(rel["af"], 1)))))
    return dict(matrix=matrix, chain_verified=led.verify_chain(),
                n_records=len(led.records),
                tcb_inventory=dict(components=["manifest_validator", "egress_monitor",
                                               "hash_chain_ledger"],
                                   loc=0, note="plaintext release-object layer; "
                                   "the sole-egress binding is built separately "
                                   "in code/enforce.py"))


if __name__ == "__main__":
    import pprint
    pprint.pprint(violation_injection_suite())
