import structure
import svg_painter
import math
import svgwrite
import os
import sys

def Draw(struct, output_file):
    root, ext = os.path.splitext(output_file)
    for idx in range(len(struct.tetrads)):
        quadruplex = structure.Quadruplex(struct, idx)
        config = svg_painter.Config(1.0) # 1.0 - Scale
        svg_maker = svg_painter.SvgMaker(config, root + "_" + str(idx) + ext, \
                quadruplex)

        # OPTIMIZE, Takes argument "optimizer" with location to the optimizer
        # binary. Default is "./svg_optimizer"
        quadruplex.Optimize()

        # Prepare + Draw
        svg_maker.DrawAll()

        # Save
        svg_maker.svg.save(pretty=True)

        # Also draw sinle tetrads from quadruplex structure.
        if (len(struct.single_tetrads[idx]) > 1):
            for tetrad_idx in range(len(struct.single_tetrads[idx])):
                quadruplex = structure.Quadruplex(struct, idx, tetrad_idx)
                config = svg_painter.Config(1.0) # 1.0 - Scale
                svg_maker = svg_painter.SvgMaker(config, root + "_" + str(idx) + \
                        "_" + str(tetrad_idx) + ext, quadruplex)

                # OPTIMIZE, Takes argument "optimizer" with location to the optimizer
                # binary. Default is "./svg_optimizer"
                quadruplex.Optimize()

                # Prepare + Draw
                svg_maker.DrawAll()

                # Save
                svg_maker.svg.save(pretty=True)


def DrawFromString(json, output_file):
    Draw(structure.Structure().fromString(json), output_file)


def DrawFromFile(filename_json, output_file):
    Draw(structure.Structure().fromFile(filename_json), output_file)

if len(sys.argv) != 3:
    print("Please provide 2 arguments. ./main.py <input_json> <output_svg>")
    exit()

DrawFromFile(sys.argv[1], sys.argv[2])
