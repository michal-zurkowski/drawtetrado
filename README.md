# Project description

DrawTetrado is a Python application to visualize quadruplexes and
G4-helices in DNA and RNA structures. It generates publication-quality
SVG files containing layer diagrams. They show the tetrads as a stack,
with each position having four nucleobases colored according to anti or
syn conformation. In addition, DrawTetrado visualizes the strands with
arrows for an accessible overview of their directionality and visual
determination of loop types.

DrawTetrado automatically optimizes the layout. The result is a visually
pleasing and readable image, even for the most complex cases like
V-loops and G4-helices (dimers).

# Installation

    pip install drawtetrado

# Dependencies

The project is written in Python 3.6+ and requires
[pycairo](https://github.com/pygobject/pycairo) and
[svgwrite](https://github.com/mozman/svgwrite). The internal
optimization routine is written in C++ and requires
[Cython](https://cython.org/) and a C++20-compliant compiler (e.g. GCC
8+).

DrawTetrado parses the output of
[ElTetrado](https://github.com/tzok/eltetrado) (Zok *et al.*, 2022;
Popenda *et al.*, 2020; Zok *et al.*, 2020). If you do not have that
installed, please install DSSR (Lu *et al.*, 2015) and run:

    pip install eltetrado

# Usage

    usage: drawtetrado [-h] [--angle ANGLE] [--tetrad-spacing TETRAD_SPACING]
                       [--nucl-longer NUCL_LONGER] [--nucl-shorter NUCL_SHORTER]
                       [--nucl-spacing NUCL_SPACING] [--line-width LINE_WIDTH]
                       [--no-label-tilted] [--no-label-nucleotide-full]
                       [--no-label-nucleotide] [--no-label-chain] [--no-label-number]
                       [--label-font LABEL_FONT] [--label-size LABEL_SIZE]
                       [--se-label-size SE_LABEL_SIZE]
                       [--se-label-spacing SE_LABEL_SPACING]
                       [--color-connection COLOR_CONNECTION]
                       [--color-border COLOR_BORDER] [--color-text COLOR_TEXT]
                       [--color-gb-anti COLOR_GB_ANTI] [--color-gb-syn COLOR_GB_SYN]
                       [--color-gb-unknown COLOR_GB_UNKNOWN]
                       input output

    positional arguments:
      input                 path to input JSON generated by ElTetrado
      output                path to output SVG file template

    optional arguments:
      -h, --help            show this help message and exit
      --angle ANGLE         Angle of the resulted drawing in degrees. Angle between
                            front and left edges of the drawing. Reasonable values are
                            probably within 70 and 40 deg. [default=50.0]
      --tetrad-spacing TETRAD_SPACING
                            Spacing between tetrad layers. Spacing of 0 will result in
                            back edge of lower layer aligning with front edge of the
                            higher layer. [default=50.0]
      --nucl-longer NUCL_LONGER
                            Length of the "longer" edge of the nucleotide box. This
                            value should be larger than --nucl-shorter to preserve
                            proper visualization of ONZ positive and negative
                            classification of the tetrad. [default=100.0]
      --nucl-shorter NUCL_SHORTER
                            Lenght of the "shorter" edge of the nucleotide box. This
                            value should be smaller than --nucl-longer to preserve
                            proper visualization of ONZ positive and negative
                            classification of the tetrad. [default=70.0]
      --nucl-spacing NUCL_SPACING
                            Spacing between individual nucleotides in the same tetrad
                            layer. [default=10.0]
      --line-width LINE_WIDTH
                            Width of the lines on the drawing. Exact witdth of
                            connection lines, tetrad border line. Points at the corners
                            of tetrads are using 70% of this value as width.
                            [default=3.0]
      --no-label-tilted     By default labels are skewed and put on labels. This means
                            the text skew is affected by the --angle of the drawing and
                            can be less readable with some configurations. Disabling it
                            will at most rotate the text instead of skewing it.
      --no-label-nucleotide-full
                            Exclude full name from the nucleotide label. Includes naming
                            like DG instead of G, BRU instead of u. Full labels can
                            provide more information about nucleotides.
      --no-label-nucleotide
                            Exclude name from the nucleotide label like A, C, T, G
                            Setting this flag also sets --no-label-nucleotide-full.
      --no-label-chain      Exclude chain information from the nucleotide label.
      --no-label-number     Exclude nucleotide index from the nucleotide label.
      --label-font LABEL_FONT
                            Font family used for all labels.
                            [default="Arial, Helvetica"]
      --label-size LABEL_SIZE
                            Font size of the nucleotide labels only. [default=20.0]
      --se-label-size SE_LABEL_SIZE
                            Font size of 5' and 3' labels only. [default=24.0]
      --se-label-spacing SE_LABEL_SPACING
                            Spacing between 5', 3' labels and nucleotide. Can put
                            further away or closer to prevent potential overlap with
                            strand between nucleotides. [default=20.0]
      --color-connection COLOR_CONNECTION
                            Color of connecting line between nucleotides in hex.
                            [default="#000000"]
      --color-border COLOR_BORDER
                            Color of rectangular border of tetrad layer in hex.
                            [default="#E23D28"]
      --color-text COLOR_TEXT
                            Color of nucleotide labels and 5' 3' labels in hex.
                            [default="#000000"]
      --color-gb-anti COLOR_GB_ANTI
                            Color of nucleotide box with anti glycosidic bond in hex.
                            [default="#FFCC38"]
      --color-gb-syn COLOR_GB_SYN
                            Color of nucleotide box with anti glycosidic bond in hex.
                            [default="#FF992B"]
      --color-gb-unknown COLOR_GB_UNKNOWN
                            Color of nucleotide box with unknown glycosidic bond in hex.
                            [default="#BB9977"]


    The output path is a template. For example, if output=/tmp/out.svg, the
    resulting files will be /tmp/out_0.svg, /tmp/out_1.svg and so on (for as many
    distinct quadruplexes as there are in the input file)

# Visual customization

DrawTetrado allows for wide veriaty of changes to the visual representation of
the resulting drawing. The presented drawing should provide a visual aid with
which arguments impact which exect aspects of the drawing.

Nucleotide labels are made from 3 parts:
```
  A.DG12
  A  - Chain of the nucleotide. If --no-label-chain is set, this information (and .)
       is not included in the creation of the labels.
  DG - Full name of the nucleotide. Setting --no-label-nucleotide-full will result
       in normal short names like A, C, T, G. --no-label-nucleotide will remove
       nucleotide name part altogether.
  12 - Index of the nucleotide. Can be removed with --no-label-number
```

Examples of how A.DG12 label would look like with different label options:
```
A.DG12 - Default
A.G12  - --no-label-nucleotide-full
A.12   - --no-label-nucleotide
DG12   - --no-label-chain
A.DG   - --no-label-number
A      - --no-label-number --no-label-nucleotide

```

![Visual changes](https://github.com/michal-zurkowski/drawtetrado/blob/main/2hy9_visuals.svg?raw=true)

# Examples

## Human telomere DNA quadruplex

![[2HY9: Human telomere DNA quadruplex structure in K+ solution hybrid-1
form](https://www.rcsb.org/structure/2hy9)](https://github.com/michal-zurkowski/drawtetrado/blob/main/2hy9.svg?raw=true)

## V-loop

![[6TCG: 2’-F-riboguanosine and 2’-F-arabinoguanosine modified
G-quadruplex with V-loop and all-syn
G-tract](https://www.rcsb.org/structure/6tcg)](https://github.com/michal-zurkowski/drawtetrado/blob/main/6tcg.svg?raw=true)

## G4-helix (dimer)

![[1MYQ: An intramolecular quadruplex of (GGA)(4) triplet repeat DNA
with a G:G:G:G tetrad and a G(:A):G(:A):G(:A):G heptad, and its dimeric
interaction](https://www.rcsb.org/structure/1myq)](https://github.com/michal-zurkowski/drawtetrado/blob/main/1myq.svg?raw=true)

# Bibliography

<div id="refs" class="references csl-bib-body">

1.  ONQUADRO: a database of experimentally determined quadruplex
    structures. T. Zok, N. Kraszewska, J. Miskiewicz, P. Pielacinska, M.
    Zurkowski, M. Szachniuk. *Nucleic Acids Research*. 2022.
    50(D1):D253–D258.
    doi:[10.1093/nar/gkab1118](https://doi.org/10.1093/nar/gkab1118)

2.  Topology-based classification of tetrads and quadruplex
    structures. M. Popenda, J. Miskiewicz, J. Sarzynska, T. Zok, M.
    Szachniuk. *Bioinformatics*. 2020. 36(4):1129–1134.
    doi:[10.1093/bioinformatics/btz738](https://doi.org/10.1093/bioinformatics/btz738)

3.  ElTetrado: a tool for identification and classification of tetrads
    and quadruplexes. T. Zok, M. Popenda, M. Szachniuk. *BMC
    Bioinformatics*. 2020. 21(1):40.
    doi:[10.1186/s12859-020-3385-1](https://doi.org/10.1186/s12859-020-3385-1)

4.  DSSR: an integrated software tool for dissecting the spatial
    structure of RNA. X.-J. Lu, H.J. Bussemaker, W.K. Olson. *Nucleic
    Acids Research*. 2015. 43(21):e142.
    doi:[f73r8c](https://doi.org/f73r8c)

</div>
