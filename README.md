# ANNOY-C
*A C implementation of Spotify’s “Approximate Nearest Neighbor Oh Yeah” (Annoy) algorithm, with plans for Python bindings.*

---

## Features

| Status | Feature |
|--------|---------|
| OK | Build-time KD-tree-style index with random hyperplane splits |
| OK | Cosine-similarity search (`searchTopK`) |
| OK | Binary serialization / deserialization of the tree |
| WIP | Python C-extension wrapper (planned) |
| WIP | SQLite metadata store (`uuid → JSON`) |
| WIP | Additional distance metrics (Euclidean, dot-product, …) |


