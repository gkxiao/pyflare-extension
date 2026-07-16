#!/usr/bin/env python
"""
decimer-segmentation.py - 使用 DECIMER 从 PDF 中分割化学结构式。

Usage:
    python decimer-segmentation.py -i input.pdf -o compounds_screenshots
    python decimer-segmentation.py -i input.pdf -o output --start-page 40 --start-compound 101 --dpi 300

输出文件命名: page-{page}-compound-{id}.png
"""

import argparse
import os
import sys
import time

os.environ["TF_ENABLE_ONEDNN_OPTS"] = "0"
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"

import fitz  # pymupdf
import numpy as np
import cv2
from decimer_segmentation import segment_chemical_structures


def main():
    parser = argparse.ArgumentParser(
        description="DECIMER-Segmentation: extract chemical structure images from PDF"
    )
    parser.add_argument("-i", "--input", required=True, help="Input PDF file path")
    parser.add_argument("-o", "--output-dir", required=True, help="Output directory for PNG screenshots")
    parser.add_argument("--start-page", type=int, default=1,
                        help="Document page number of the first PDF page (default: 1)")
    parser.add_argument("--start-compound", type=int, default=1,
                        help="Starting compound ID number (default: 1)")
    parser.add_argument("--dpi", type=int, default=300,
                        help="Render DPI for PDF pages (default: 300)")
    parser.add_argument("--prefix", type=str, default="",
                        help="Optional prefix for output filenames")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"ERROR: input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    doc = fitz.open(args.input)
    n_pages = doc.page_count
    print(f"Input: {args.input}")
    print(f"Pages: {n_pages}")
    print(f"Output: {args.output_dir}")
    print(f"DPI: {args.dpi}")
    print(f"Start page: {args.start_page}, start compound: {args.start_compound}")
    print(f"{'='*60}")

    page_num = args.start_page
    compound_id = args.start_compound
    total_saved = 0
    t_total = time.time()

    for i in range(n_pages):
        t0 = time.time()

        page = doc[i]
        mat = fitz.Matrix(args.dpi / 72, args.dpi / 72)
        pix = page.get_pixmap(matrix=mat)
        img = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)

        if img.shape[2] == 4:
            img = cv2.cvtColor(img, cv2.COLOR_RGBA2BGR)
        elif img.shape[2] == 1:
            img = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)

        segments = segment_chemical_structures(img, expand=True)
        dt = time.time() - t0

        print(f"Page {page_num}: {len(segments)} structures ({dt:.1f}s)")

        for seg in segments:
            prefix = f"{args.prefix}" if args.prefix else ""
            fname = f"{prefix}page-{page_num}-compound-{compound_id}.png"
            fpath = os.path.join(args.output_dir, fname)

            if len(seg.shape) == 2:
                seg = cv2.cvtColor(seg, cv2.COLOR_GRAY2BGR)
            elif seg.shape[2] == 4:
                seg = cv2.cvtColor(seg, cv2.COLOR_RGBA2BGR)

            cv2.imwrite(fpath, seg)
            compound_id += 1
            total_saved += 1

        page_num += 1

    doc.close()
    elapsed = time.time() - t_total
    print(f"\n{'='*60}")
    print(f"Done. {total_saved} structures in {elapsed:.0f}s ({elapsed/n_pages:.1f}s/page)")
    print(f"Output: {args.output_dir}")


if __name__ == "__main__":
    main()
