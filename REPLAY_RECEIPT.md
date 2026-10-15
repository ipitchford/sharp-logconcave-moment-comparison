# Clean replay receipt

Both actual local runs completed all eleven production certificates and semantic output checks. This is same-code replay, not independent reimplementation or formal verification.

| Mode | Certificates passed | Sum of per-command seconds | Portable evidence |
|---|---:|---:|---|
| normal | 11 | 117.80 | [receipt](receipts/replay-normal.json) |
| optimized | 11 | 117.47 | [receipt](receipts/replay-optimized.json) |

The eleven conventional certificate files are the regenerated normal-mode outputs, byte-for-byte. Optimized outputs and raw logs remain privately retained; their hashes are recorded. Each receipt binds the replayed code and tests, including the shared endpoint serializer. New packaging/CI files were added after these local scientific runs and are not claimed as replay inputs.

Portable replay: `bash code/reproduce.sh /path/to/python /path/to/NEW_OUTPUT [--optimized]`. Install `requirements.txt` in an isolated environment first. Run `python -B checks/test_certificates.py /path/to/NEW_OUTPUT` (also `-O` for optimized tests). Use a new output directory outside the checkout.

Both six-method regression suites passed with actual certificates. Eight deliberately false bulk-control runs and eighteen malformed-input/overwrite requests were rejected. See [semantic-control receipt](receipts/semantic-controls.json) for actual result and private-log hash bindings.

CI status must be read from the actual GitHub run; this receipt makes no Linux execution or public publication claim.
