# Proyecto Series de Tiempo
## Concentración Ambiental de Contaminantes Atmosféricos en La Florida

Fecha: 22/10/2026
Integrantes: Tomás Pérez, Julián Vargas, Benjamin R, Luna Garcés, Pilar Vallotton

### Contenido:
Código fuente y base de datos del proyecto de series de tiempo relacionado a la concentración ambiental de contaminantes atmosféricos en la comuna de la florida.

### Propuesta del Proyecto:
Se busca realizar un estudio medioambiental con el uso de series de tiempo e información de la concentración ambiental de contaminantes atmosféricos de la estación ubicada en la comuna de La Florida. El compuesto a estudiar será el Material Particulado Fino MP2,5 el cual se encuentra representado como promedios móviles de 24 horas, basados en el monitoreo continuo. Este compuesto se asocia a los niveles de concentración que originan situaciones de emergencia ambiental, para aquellas estaciones que cuentan con monitoreo en línea.
El Material Particulado Fino (MP2,5) es un contaminante atmosférico que corresponde a partículas líquidas o sólidas en suspensión con un diámetro aerodinámico menor a 2,5μm, el cual, debido a su tamaño, puede recorrer todo el tracto respiratorio e incluso atravesar los alvéolos pulmonares, lo que presenta un alto riesgo para la salud de la población (más dañino para la salud que el MP10). Su composición química es variable y dependerá de sus fuentes de emisión.
Las principales fuentes de emisión del MP2,5 son variadas y pueden ser de origen antropogénico directo (quema de combustibles, refinería de metales, incendios forestales, etc.), indirecto (aerosoles, reacciones químicas, etc.) y natural.
Este tipo de contaminantes atmosféricos provoca efectos adversos en el ser humano como aumento en la mortalidad prematura, aumento en admisiones hospitalarias en pacientes con enfermedades pulmonar obstructiva crónica y con enfermedad cardiovascular, exacerbación de los síntomas e incremento del asma, aumento de riesgo de infartos al miocardio, inflamación pulmonar, inflamación sistémica, disfunciones endoteliales y vasculares, desarrollo de aterosclerosis, incremento en la incidencia de infecciones y cáncer respiratorio.

### Fuente de los Datos:
Sistema de Información Nacional de Calidad del Aire (SINCA), Ministerio del Medio Ambiente (https://sinca.mma.gob.cl), Región Metropolitana, estación La Florida, parámetro MP2,5, registros diarios validados. ([Gráfico MP2,5 La Florida](https://sinca.mma.gob.cl/cgi-bin/APUB-MMA/apub.htmlindico2.cgi?page=pageFrame&header=La%20Florida%20(Acreditada)&macropath=./RM/D12/Cal/PM25&macro=PM25.diario.diario&from=020101&to=260823&))

### Series a analizar:
Concentración promedio mensual de material particulado fino (MP2,5), en μg/m³ dentro de la comuna de La Florida

### Frecuencia y periodo:
Datos mensuales entre enero de 2002 y agosto de 2026, con 295 observaciones.

### Pregunta a responder:
¿Existen patrones de tendencia y estacionalidad en la concentración promedio mensual de MP2,5 en La Florida, y qué tan bien pueden predecirse sus valores mediante modelos de series de tiempo?
