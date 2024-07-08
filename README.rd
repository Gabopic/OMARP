OMARP (OMA reduction pipeline)

1)Introducción:
-MAS500 es un telescopio de 50cm de apertura con 6500 mm de distancia focal, con un FoV de 24.57'. Se ubica en el observatorio el Sauce, comuna Río Hurtado, Cuarta región, Chile a 1600 metros sobre el nivel del mar.
-Lista de filtros (Sloan:(u,g,i,r,z) y Jhonson/Bessel:(V y B))

2)Software: 
-Instalación:Usar pip
-ESte programa cumple la función de reducir los datos de las imágenes obetindas de MAS500 para limpiarlas y además hacer las calibraciones necesarias (Astrometría y Fotometría).

3)Input:

-BPM: Para la creación del bad pixel map se usaron darks de destinos tiempos de exposición para estudiar la distribución de ADUs en cada uno, los tiempos fueron de 50, 100, 200, 300 segundos tomados el 6 de Junio de 2024, 
se calculó el dark current en cada grupo de imágenes y se compararon los resultados de cada tiempo de exposición para ver qué puntos se comportan de manera lineal.

-Campos de flat (tbw)

-Raw,Bias,Darks,Flats

4)Uso:
-Para hacer las reducciones hay dos funciones principales, tables() y Redux()
-Mostrar fórmulas

5)Consideraciones:
-Astrometría (Mencionar líneas faltantes en el header)
-Fotometría (zeropoint) (En desarrollo)
-Estabilidad de los darks

6)Output: 
-Fits
-Imágenes
-Ejemplo
