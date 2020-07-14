# Bienvenidos a metodos RK(s) y RKN(s) en MatLab!

Este repositorio es el resultado de un trabajo de investigación para obtener el titulo de Magister en Matemáticas. Aquí encontraras implementación de los métodos numéricos  **Runge-Kutta** y  **Runge-Kutta-Nystrom** tanto en su presentación explicita como implícita. 

## Metodos implementados

Los métodos implementados que encontraras en este repositorio son los siguientes: 

 1. Metodos Runge-Kutta: 
	 - Explicitos:
		 - Euler (RK 1er orden)
		 - Euler mejorado (RK 2do orden)
		 - Heun (RK 3er orden)
		 - Runge-Kutta (RK 4to orden)
	- Implicitos:
		- Gauss 1 etapa (RK 2do orden)
		- Gauss 2 etapas (RK 4to orden)
		- Gauss 3 etapas (RK 6to orden)
		- Kuntemann y Butcher (RK 8vo orden)
 2. Metodos Runge-Kutta-Nystrom:
	 - Explicitos:
		 - RKN:
			- Directo (RKN 3er orden)
			- Transformado (RKN 3er orden)
			- Nystrom (RKN 4to orden)
		 - RKN simple:
			 - Nystrom simple (RKN 4to orden)
			 - Nystrom simple orden 5 (RKN 5to orden)
			 - Nystrom especial orden 5 (RKN 5to orden)
			 - DOPRI5(4)7FM (RKN 7mo orden)
	 - Implicito:
		 - Runge-Kutta-Nystrom implicito 8vo orden

**Nota:** Estos metodos se encuentran repetidos debido a su jerarquia en el repositorio.

## Jerarquia de carpetas y subcarpetas

La jerarquia está organizada en primer lugar por el tipo de metodo, luego por el tipo de ejercicio y por último el tipo de presentación (explicito o implcito). Al ubicarse en el tipo de presentación ya se encuentra el codigo con una solución particular.

Foto

Los tipos de ejercicios son:

 1. Ecuacion diferencial ordinario de primer orden
 2. Sistema de ecuaciones diferenciales ordinaria de primer orden
 3. Ecuacion diferencial ordinaria de segundo orden $f(x,y)$
 4. Sistema de ecuaciones diferenciales ordinaria de segundo orden $f(t,x,y)$ y $g(t,x,y)$

## Ejercicios resueltos

Está implementación de cada metodo cuenta con la solución de un ejercicio en especifico, pero el codigo es totalmente gratis y puede ser descargado para posteriormente ser modificado y usado con cualquier tipo de finalidad propia.

 - Ecuacion diferencial ordinaria de primer orden
$$
\begin{cases}
 y'=\cos(\sqrt{y-x}) \\
 y(0)=\pi^{2}
 \end{cases}
 $$
 
 - Sistema de ecuaciones diferenciales ordinaria de primer orden
$$
 \begin{cases}
 x'(t) = \frac{\cot{(t)}}{2e^{\sqrt{t}} \tan{(t^2)}}  + \frac{1}{2t}x(t) \\
 y'(t) = \frac{2t}{e^{\sqrt{t}}} y^2(t) +\frac{\sqrt{\sin{(t)}}}{2}\frac{y(t)}{x(t)}+2\csc{(t)} e^{\sqrt{t}}x^2(t) \\
 x(0.01)=0.0099999166667361106977490148082 \\
 y(0.01)=0.0001105170921759550699799478876
 \end{cases}
 $$
- Ecuacion diferencial ordinaria de segundo orden $f(x,y)$
$$
 \begin{cases}
 y'' = \sin(x-y^2)  \\
 y'(0) = 1 \\
 y(0) = 0
 \end{cases}
$$
 - Sistema de ecuaciones diferenciales ordinaria de segundo orden $f(t,x,y)$ y $g(t,x,y)$
$$
\begin{cases}
 x'' = \frac{2y}{x^2 + y^2} - 4t^2x  \\
 y'' = -2x - \frac{4t^2y}{x^2 + y^2} \\
 x(0) = 0 \hspace{3mm} x'(0) = 0 \\
 y(0) = 1 \hspace{3mm} y'(0) = 0
 \end{cases}
$$

## Tesis