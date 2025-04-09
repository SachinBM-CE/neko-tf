#!/usr/bin/env python3
"""
Python Script: Convert a Set of PNG Images to an MP4 Video

This script:
1. Searches for all PNG images in the specified directory.
2. Sorts the images in alphabetical (or numeric) order.
3. Creates a video clip from the image sequence at a specified frame rate.
4. Exports the video as an MP4 file using the 'libx264' codec.

Usage:
    python png_to_mp4.py

Requirements:
    - Python 3.x
    - moviepy (install via pip: pip install moviepy)
"""

import os
import glob
from moviepy import ImageSequenceClip

def create_video_from_pngs(image_folder=".", output_file="output.mp4", fps=24):
    """
    Create an MP4 video from PNG images in a folder.

    Parameters:
        image_folder (str): Path to the folder containing PNG images.
        output_file (str): Filename for the output video.
        fps (int): Frames per second for the output video.
    """
    # Construct the search path and sort the files
    search_path = os.path.join(image_folder, "*.png")
    image_files = sorted(glob.glob(search_path))

    # Check if any images were found
    if not image_files:
        print(f"No PNG images found in the folder: {image_folder}")
        return

    print(f"Found {len(image_files)} image(s). Converting to video...")

    # Create a video clip from the image sequence
    clip = ImageSequenceClip(image_files, fps=fps)

    # Write the result to a file
    clip.write_videofile(output_file, codec='libx264')
    print(f"Video saved as {output_file}")

if __name__ == "__main__":
    # Customize parameters if needed:
    IMAGE_FOLDER = "/tmp/sachinbm/neko-tf/neko/examples/0_RESULTS/maggradu/"   # Current directory; change if images are elsewhere
    OUTPUT_FILE = "/tmp/sachinbm/neko-tf/neko/examples/0_RESULTS/maggradu/output.mp4"
    FPS = 24

    create_video_from_pngs(image_folder=IMAGE_FOLDER, output_file=OUTPUT_FILE, fps=FPS)
