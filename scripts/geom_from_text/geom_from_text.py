import cadquery as cq

def create_text_cutout_brep(text: str, filename: str = "text_cutout.brep"):
    """
    Erstellt eine rechteckige Platte, schneidet den übergebenen Text aus
    und speichert das Ergebnis als BREP-Datei.
    """
    # Parameter
    
 #   plate_width = 200.0   # Breite des Rechtecks (X-Achse)
 #    plate_height = 80.0   # Höhe des Rechtecks (Y-Achse)
    plate_thickness = 30.0 # Dicke der Platte (Z-Achse)
    padding = 20 # Dicker der Platte um den Text rum    

    # Parameter für den Text
    font_size = 30.0      # Schriftgröße
    font_thickness = 3*plate_thickness # Extrusionstiefe des Textes (muss dicker als die Platte sein für sauberen Durchbruch)


    # Erstelle den Text als 3D-Objekt
    # Wir zentrieren den Text standardmäßig auf der XY-Ebene
    text_geometry = (
        cq.Workplane("XY")
        .text(text, font_size, font_thickness, font="Arial", halign="center", valign="center")
    )
    # Verschiebe den Text in Z Achse, damit er die Platte durchschneidet
    text_geometry = text_geometry.translate ((0,0,-plate_thickness))

    # Erstelle die Basis-Grundplatte (Rechteck)
    bounding_box = text_geometry.val().BoundingBox()
    text_width = bounding_box.xlen
    text_height = bounding_box.ylen
    plate_width = text_width + (2 * padding)
    plate_height = text_height + (2 * padding)

    plate = cq.Workplane("XY").box(plate_width, plate_height, plate_thickness)
    
    # Schneide den Text aus der Platte aus (Boolesche Differenz)
    result = plate.cut(text_geometry)

    # Exportiere das Ergebnis als BREP-Datei
    cq.exporters.export(result, filename, cq.exporters.ExportTypes.BREP)
    print(f"Erfolgreich exportiert: {filename}")

# Beispielaufruf
if __name__ == "__main__":
    mein_text = "Sehr cool! :-)"
    print ("Creating brep file for text: \"" + mein_text + "\"")
    create_text_cutout_brep(mein_text, "text.brep")
