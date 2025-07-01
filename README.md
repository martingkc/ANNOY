# ANNOY-C Vector Index

A minimal C implementation of [ANNOY (Approximate Nearest Neighbors Oh Yeah)](https://github.com/spotify/annoy), the high-dimensional index behind Spotify’s music recommendations.  
This version builds a binary, on-disk tree of 512-dimensional vectors, supports fast cosine-similarity search, and can save/load the index to/from a file.

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

## References 
-  [Erik Bernhardsson, “Nearest neighbors and vector models, part 2”](https://erikbern.com/2015/10/01/nearest-neighbors-and-vector-models-part-2-how-to-search-in-high-dimensional-spaces.html)
-  [ANNOY (Approximate Nearest Neighbors Oh Yeah)](https://github.com/spotify/annoy)
