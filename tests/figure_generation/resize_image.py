from PIL import Image
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent

IMAGES_DIR_PATHS = SCRIPT_DIR / 'original'
OUTPUT_DIR_PATHS = SCRIPT_DIR / 'resized'


def resize_image(image_path):
    # Load your image
    img = Image.open(image_path)

    # Original size
    # original_size = img.size  # (width, height)

    # # Scale factor — e.g., 2x resolution
    # scale_factor = 6
    # new_size = (original_size[0] * scale_factor, original_size[1] * scale_factor)
    #
    # # Resize image
    # resized_img = img.resize(new_size, Image.LANCZOS)

    # Save with higher DPI (e.g., 300)
    img.save(OUTPUT_DIR_PATHS / image_path.name, dpi=(600, 600))
    print(f"Resized {image_path.name} and saved to {OUTPUT_DIR_PATHS / image_path.name}")


def main():
    # Create output directory if it doesn't exist
    OUTPUT_DIR_PATHS.mkdir(parents=True, exist_ok=True)

    # Iterate over all images in the directory
    for image_path in IMAGES_DIR_PATHS.glob('*.png'):
        resize_image(image_path)


if __name__ == "__main__":
    main()
