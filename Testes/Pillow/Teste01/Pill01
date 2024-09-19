from PIL import Image
import os

def get_dominant_color(image_path):
    with Image.open(image_path) as img:
        img = img.convert("RGB")
        colors = img.getcolors(img.size[0] * img.size[1])
        # Retorna a cor mais frequente
        dominant_color = max(colors, key=lambda item: item[0])[1]
        return dominant_color

def create_folder_if_not_exists(folder_path):
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

def move_image_by_color(image_path, output_folder, dominant_color):
    r, g, b = dominant_color
    # Converte a cor para um formato simples para usar no nome da pasta
    color_folder_name = f"{r}_{g}_{b}"
    color_folder_path = os.path.join(output_folder, color_folder_name)
    create_folder_if_not_exists(color_folder_path)
    
    # Move a imagem para a pasta correspondente
    image_name = os.path.basename(image_path)
    new_image_path = os.path.join(color_folder_path, image_name)
    os.rename(image_path, new_image_path)

# Caminho para as imagens e a pasta de saída
input_folder = "path/to/your/images"
output_folder = "path/to/separated/images"

for filename in os.listdir(input_folder):
    if filename.endswith((".jpg", ".jpeg", ".png")):
        image_path = os.path.join(input_folder, filename)
        dominant_color = get_dominant_color(image_path)
        move_image_by_color(image_path, output_folder, dominant_color)
