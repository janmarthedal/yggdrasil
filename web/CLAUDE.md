# Posts and media on the Finite Element Method and the Yggdrasil project

## Purpose
The folder ./posts/ contains a series of posts introducing the finite element
method, accompanying the `yggdrasil` library in this repository.
The posts explain the mathematical foundation of the finite element method
and the code constructs that make up the library in relation to the theory.
There will be little abstract mathematical theory, but references to the
underlying theory should be included. Resources accessible online are preferred,
but authoritative books and papers may also be mentioned.

**Each post must stand on its own as finite element theory.** A reader who never
opens the `yggdrasil` source — or who skips every library reference — should still
be able to follow the post in full. Library references *exemplify* the theory in
concrete code; they are never load-bearing for understanding it.

It is an important goal that every function and class from the library are
referenced from at least one post at some point. This goal is compatible with the
principle above: the reference shows the theory realized in code, but the
surrounding prose carries the explanation.

## Post guidelines
- Plenty of examples, illustrations and animations should be included.
- Each post is of small size and as self-contained as possible.
- Use a shared set of notation, which is listed in `notation.md`.
  Update this file as posts are written.

## Writing conventions
- Use markdown syntax.
- Reference library source files using links, e.g., reference
  yggdrasil/assemble.py as "[`assemble.py`](LIBROOT/assemble.py)".

## Library references
- Theory first: write the post so that a reader who skips every library mention
  still understands it completely. Develop each concept in plain mathematical
  prose before pointing to the code that implements it.
- Keep library references *inline but skippable*: place them near the relevant
  theory, but confine them to their own self-contained sentence(s) or paragraph(s)
  that can be deleted without breaking the surrounding text. Do not thread class
  names, method names, or array-shape details through the sentences that carry the
  mathematical argument.
- Frame each reference as an example, not a definition. Prefer connective openings
  like "In the library, this is …" or "Yggdrasil realizes this as …" so it reads
  as a concrete illustration of the theory just presented.
- Do not introduce headings for these asides (the no-headings convention still
  applies); use the connective framing above instead.
- Reference image files by path relative to `MEDIAROOT`.
- Reference other posts by path relative to `POSTROOT`. The post link for file
  `some-post.md` should be `POSTROOT/some-post/`.
- Math is written in LaTeX delimited by `$...$` (inline) and `$$...$$` (display).
- Prefer SVG over PNG for inline images when the illustration is fit for vector
  graphics.
- Consider using [three.js](https://threejs.org/) for 3D illustrations.
- Headings should be avoided. Let the text flow as continuous prose.

## Post Specification
The file [post specification](./post-specification.md) contains a list
of current and planned posts in chronological order.
In addition to each post's filename and title, it also includes a specification
of each post's content.
