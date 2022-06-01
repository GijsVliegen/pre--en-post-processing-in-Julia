Deze repository dient voor de wetenschappelijke vorming van Tim Olde en Gijs vliegen, met als promotor Martin Diehl
project: Pre- and postprocessing in Julia

Structuur van deze repository

* Julia
  * rotations.jl 			De Julia-implementatie van damask.Rotation.py
  * rotationsTest.jl			Test-suite
  * rotations_timing.jl		Experimenten om de efficiëncy te meten van de Julia-implementatie
* Python
  * rotations_DAMASK_timing.py	Experimenten om de efficiëncy te meten van de DAMASK-implementatie
  * rotations_naïef.py			De Naïeve implementatie van damask.Rotation.py 
  * rotations_naïef_timing.py   	Experimenten om de efficiëncy te meten van de naïeve implementatie
* .gitignore
* README.md
* charCount.py				Script om de leesbaarheid te meten van een stuk code
