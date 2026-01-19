# Implementation of Hierarchical Z-Buffer Algorithm
Computer Graphics Course Assignment (Zhejiang University, Fall/Winter Semester 2025)

## Features
- **Scanline:** A specially designed hash-based skip list with contiguous memory layout, supporting low-constant O(1) insert/delete and O(n) active-list traversal.

  *10x faster than `std::multimap`.*

- **Quadtree Hierarchy:** A lazy-update hybrid strategy. hierarchical at the patch level, with scanline rendering within each patch. The quadtree is updated only after a certain number of accumulated Z-buffer replacement failures.

  *4× speedup on small scenes, approaching the performance of standard scanline. 5x faster than scanline on large scene.*

- **Octree Hierarchy:** A counting-sort–accelerated, lazily built octree. Counting sort reduces memory accessesand keeps each subtree contiguous in memory. Child nodes are created only when the traversal reaches the node during rendering.

  *2x faster than temporary `std::vector` + prebuilt octree. 7x faster than scanline on large scene.*

## Experiments
| SCENE | FACES | NAIVE | SCANLINE | QUADTREE | OCTREE |
|------|------|------|--------|--------|--------|
| Veach-Mis | 2,218 | **1.21ms** | 1.38ms | 6.98ms | 6.14ms |
| Cornell-Box | 26,678 | 7.99ms | **6.33ms** | 15.37ms | 14.46ms |
| Living-Room | 143,175 | **24.14ms** | 31.65ms | 58.11ms | 52.27ms |
| External-Wall | 566,131 | **56.69ms** | 163.31ms | 173.49ms | 123.19ms |
| Random-Tris | 1,000,000 | 330.01ms | 720.96ms | 224.63ms | **181.45ms** |
| Overlap-Planes | 2,112,500 | 1275.76ms | 2146.33ms | 412.35ms | **304.47ms** |

**Comparison with two baselines (from Github):**

Thanks to these open-source projects for providing many helpful references.

- https://yaelcassini.github.io/2023/05/16/ZBuffer-Report/

  | SCENE      | FACES   | ALGORITHM | BASELINE (ms) | OURS (ms) |
  | ---------- | ------- | --------- | ------------- | ---------------- |
  | torus1k    | 625     | Naive     | 170           | **1.86**         |
  | torus1k    | 625     | Scanline  | 21            | **1.81**         |
  | torus1k    | 625     | Quadtree  | 235           | **9.36**         |
  | torus1k    | 625     | Octree    | 226           | **9.79**         |
  | knob4k     | 3,984   | Naive     | 196           | **6.29**         |
  | knob4k     | 3,984   | Scanline  | 47            | **5.55**         |
  | knob4k     | 3,984   | Quadtree  | 496           | **10.21**        |
  | knob4k     | 3,984   | Octree    | 352           | **12.15**        |
  | teapot15k  | 15,704  | Naive     | 148           | **5.76**         |
  | teapot15k  | 15,704  | Scanline  | 311           | **7.25**         |
  | teapot15k  | 15,704  | Quadtree  | 263           | **14.61**        |
  | teapot15k  | 15,704  | Octree    | 250           | **13.31**        |
  | spiral120k | 120,156 | Naive     | 285           | **11.73**        |
  | spiral120k | 120,156 | Scanline  | 1,456         | **21.80**        |
  | spiral120k | 120,156 | Quadtree  | 441           | **32.98**        |
  | spiral120k | 120,156 | Octree    | 365           | **26.65**        |
  | bunny144k  | 144,046 | Naive     | 134           | **12.54**        |
  | bunny144k  | 144,046 | Scanline  | 1,726         | **22.41**        |
  | bunny144k  | 144,046 | Quadtree  | 407           | **45.25**        |
  | bunny144k  | 144,046 | Octree    | 272           | **35.71**        |

- https://github.com/Wajov/HiddenSurfaceRemover

  | SCENE     | FACES   | ALGORITHM | BASELINE (ms) | OURS (ms) |
  | --------- | ------- | --------- | ------------- | ---------------- |
  | cube      | 12      | Naive     | 111           | **1.79**         |
  | cube      | 12      | Scanline  | 81            | **1.64**         |
  | cube      | 12      | Quadtree  | 647           | **8.71**         |
  | cube      | 12      | Octree    | 661           | **8.06**         |
  | suzanne   | 968     | Naive     | 82            | **1.82**         |
  | suzanne   | 968     | Scanline  | 61            | **1.59**         |
  | suzanne   | 968     | Quadtree  | 585           | **8.02**         |
  | suzanne   | 968     | Octree    | 555           | **7.71**         |
  | bunny     | 69,630  | Naive     | 256           | **10.61**        |
  | bunny     | 69,630  | Scanline  | 254           | **14.63**        |
  | bunny     | 69,630  | Quadtree  | 1,107         | **29.41**        |
  | bunny     | 69,630  | Octree    | 1,001         | **25.81**        |
  | armadillo | 212,574 | Naive     | 327           | **19.08**        |
  | armadillo | 212,574 | Scanline  | 446           | **36.48**        |
  | armadillo | 212,574 | Quadtree  | 983           | **88.53**        |
  | armadillo | 212,574 | Octree    | 1,358         | **74.24**        |

## Build & Run
Tested on Linux and Windows. No support for MacOS.

```bash
python main.py <algorithm> <test_scene>
```
- Algorithm: `naive`, `scanline`, `hierarchy`, or `octree`
- Scene: path to the scene folder

**Example:** 
```bash
python main.py octree example-scenes-cg25/living-room
```

*The output `scene.bmp` will be generated in the scene folder. No frontend is provided (AI can implement it efficiently). I hate CMake, using Python instead. On Windows, you may compile with MSVC; you need to replace `__builtin_clz` with `__lzcnt`.
