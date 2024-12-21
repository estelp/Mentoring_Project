import qrcode
import os

# Les liens pour lesquels générer les codes QR
urls = {
    "MindMap_QR": "https://www.mindmeister.com/app/map/3489792035?t=sMvnsWUbpm",
    "Github_QR": "https://github.com/estelp/Mentoring_Project/tree/main"
}

# Répertoire de sauvegarde
output_dir = '/home/name/Documents/Projet_CIBiG/Mentoring_Project/QR_Code'
os.makedirs(output_dir, exist_ok=True)  # Créer le répertoire s'il n'existe pas

# Générer et sauvegarder un QR code pour chaque URL
for name, url in urls.items():
    # Créer un objet QR Code
    qr = qrcode.QRCode(version=1, box_size=10, border=4)
    qr.add_data(url)
    qr.make(fit=True)

    # Créer l'image du code QR
    img = qr.make_image(fill='black', back_color='white')

    # Spécifier le chemin de sauvegarde
    save_path = os.path.join(output_dir, f"{name}.png")

    # Sauvegarder l'image générée
    img.save(save_path)
    print(f"Le code QR pour '{name}' a été généré et sauvegardé sous '{save_path}'.")

