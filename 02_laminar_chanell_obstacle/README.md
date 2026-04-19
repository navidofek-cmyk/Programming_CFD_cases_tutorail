# 02_laminar_chanell_obstacle

Malý edukativní C++20 solver inspirovaný strukturou jednoduchých CFD solverů:

- 2D nestlačitelné laminární Navier-Stokesovy rovnice
- cell-centered finite volume na strukturované síti
- segregovaný SIMPLE-like pseudo-time postup
- obdélníkový kanál se čtvercovou překážkou
- standard library only

Nejde o produkční solver. Cílem je čitelnost, rozdělení do rozumných modulů a snadné experimentování.

## Struktura projektu

- `config.hpp` - fyzikální a numerické parametry
- `mesh.hpp`, `mesh.cpp` - strukturovaná síť, souřadnice středů buněk, maska překážky
- `field.hpp` - jednoduchý 2D wrapper nad `std::vector<double>`
- `boundary.hpp`, `boundary.cpp` - inlet, outlet, walls, obstacle boundary conditions
- `solver.hpp`, `solver.cpp` - SIMPLE-like iterační řešič
- `io.hpp`, `io.cpp` - export polí do CSV
- `main.cpp` - sestavení case, spuštění solveru, zápis výstupů
- `Makefile` - build a run

## Fyzika a case

Geometrie:

- výška překážky `H = 1`
- doména `Lx = 32H`, `Ly = 8H`
- překážka:
  - `x in [8H, 9H]`
  - `y in [3.5H, 4.5H]`

Proudění:

- vstupní rychlost `U_in = 1`
- Reynoldsovo číslo `Re = 60`
- viskozita z:

```text
nu = U_in * H / Re
```

## Numerická myšlenka

Solver používá jednoduchý segregovaný postup:

1. z aktuálního tlaku se spočítá predikce rychlostí
2. tlaková korekce se získá z Poissonovy rovnice odvozené z divergence predikovaných rychlostí
3. tlak i rychlosti se opraví
4. celý postup se opakuje v pseudo-čase, dokud se pole neustálí

Použitá schémata:

- konvekce: first-order upwind
- difuze: central difference
- tlaková korekce: iterativní Gauss-Seidel-like sweep
- stabilita: under-relaxation
  - `alpha_u = 0.5`
  - `alpha_v = 0.5`
  - `alpha_p = 0.3`

## Build

```bash
cd /home/ivand/projects/learning_cpp/cfd/programming_cfd_cases/02_laminar_chanell_obstacle
make
```

## Run

```bash
./channel_obstacle
```

nebo:

```bash
make run
```

## Výstupy

Program zapisuje:

- `u.csv`
- `v.csv`
- `p.csv`

Během běhu vypisuje:

- iteraci
- momentum residual pro `u`
- momentum residual pro `v`
- residual tlakové korekce
- maximální velikost rychlosti

## Co čekat fyzikálně

Za čtvercovou překážkou by se měla vytvořit recirkulační oblast:

- těsně za hranou překážky dojde k odtržení proudnic
- v wake oblasti vznikne zóna zpětného nebo velmi slabého proudění
- při `Re = 60` je proudění stále laminární a stacionární
- recirkulace by měla být symetrická vůči ose kanálu jen přibližně, protože překážka není na středu výšky domény úplně “samotná”, ale v konečném kanálu s horní a dolní stěnou

## Poznámka

Tento solver je záměrně zjednodušený:

- collocated arrangement
- jednoduchá tlaková korekce
- bez turbulence
- bez přesného průtokového vyvážení na profesionální úrovni
- překážka je stair-step maska na kartézské síti
- pseudo-time SIMPLE-like iterace je edukativní aproximace, ne plný průmyslový SIMPLE solver

Je vhodný hlavně pro studium organizace CFD kódu a základní pressure-velocity coupling logiky.
