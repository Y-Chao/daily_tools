import matplotlib.font_manager as fm

# Print out the font name for matplotlib
font_names = sorted(set(f.name for f in fm.fontManager.ttflist))
print("Font Names:\n", ', '.join(font_names))
