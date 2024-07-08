Documentación OMA:

1)Introducción:
-MAS500 es un telescopio de 50cm de apertura con 
-Lista de filtros (Sloan:(u,g,i,r,z) y Jhonson/Bessel:(V y B))

2)Software: 
-Instalación:Usar pip
-Explicar para qué sirve el programa (reducción, calibración)

3)Input:

-BPM: Para la creación del bad pixel map se usaron darks de destinos tiempos de exposición para estudiar la distribución de ADUs en cada uno, los tiempos fueron de 50, 100, 200, 300 segundos tomados el 6 de Junio de 2024, 
se calculó el dark current en cada grupo de imágenes y se compararon los resultados de cada tiempo de exposición para ver qué puntos se comportan de manera lineal.

-Campos de flat (tbw)

-Raw,Bias,Darks,Flats

4)Uso:
-Explicar funciones para generar las tablas, la reducción, calibración y weight map (mosaicos)
-Mostrar fórmulas

5)Consideraciones:
-Astrometría (Mencionar líneas faltantes en el header)
-Fotometría (zeropoint) (En desarrollo)
-Estabilidad de los darks

6)Output: 
-Fits
-Imágenes
-Ejemplo
