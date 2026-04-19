# 01_cavity_case

Jednoduchý edukativní C++20 příklad pro 2D lid-driven cavity flow:

- nestlačitelné Navier-Stokesovy rovnice
- uniformní kartézská síť
- konečné diference na strukturované mřížce
- projection / pressure-correction přístup
- explicitní krok pro hybnost
- jednoduché iterační řešení Poissonovy rovnice pro tlak

Cílem tohoto příkladu je čitelnost a snadná úprava, ne produkční CFD solver.

## Obsah složky

- `main.cpp` - celý solver v jednom souboru
- `run.sh` - překlad a spuštění jedním příkazem
- `cavity` - přeložený program po kompilaci
- `u.csv` - pole rychlosti ve směru `x`
- `v.csv` - pole rychlosti ve směru `y`
- `p.csv` - tlakové pole
- `cavity.vtk` - výstup pro ParaView

## Struktura case

Prakticky je case rozdělený takto:

- `main.cpp`
  - uživatelské parametry na začátku souboru
  - malé pomocné funkce pro indexaci, okrajové podmínky, Poissonův solver a export
  - hlavní časová smyčka v `main()`
- `run.sh`
  - jednoduchý wrapper pro překlad a spuštění
- výstupy `csv` a `vtk`
  - vznikají po doběhu solveru
  - mohou být znovu přepsány při dalším spuštění

## Co úloha řeší

Jde o klasický cavity case:

- čtvercová dutina `1 x 1`
- horní stěna se pohybuje konstantní rychlostí `U_lid`
- ostatní stěny jsou nehybné
- proudění je 2D a nestlačitelné

Tento problém je velmi běžný jako první test CFD kódu, protože:

- geometrie je jednoduchá
- okrajové podmínky jsou přehledné
- vzniká uvnitř dutiny charakteristická recirkulace
- dobře se na něm učí vazba mezi rychlostí a tlakem

## Fyzikální měřítka

Pro cavity flow se často sleduje Reynoldsovo číslo:

```text
Re = U_lid * L / nu
```

kde:

- `U_lid` je rychlost horního víka
- `L` je charakteristický rozměr dutiny
- `nu` je kinematická viskozita

V tomto case je:

- `L = 1`
- výchozí `U_lid = 1`
- výchozí `nu = 0.1`

Takže výchozí hodnota je:

```text
Re = 1 * 1 / 0.1 = 10
```

To je poměrně mírné proudění, což je vhodné pro jednoduchý explicitní edukativní solver.

## Numerická metoda

Použitý postup je jednoduchá projection metoda:

1. Z aktuálních rychlostí se sestaví pravá strana tlakové Poissonovy rovnice.
2. Vyřeší se tlak z Poissonovy rovnice iterativně.
3. Rychlosti se aktualizují explicitně pomocí konvekce, difuze a tlakového gradientu.
4. Znovu se aplikují okrajové podmínky.

Použité zjednodušení:

- kolokovaná mřížka
- explicitní časový krok
- jednoduchý iterativní tlakový solver
- bez pokročilé kontroly konvergence
- bez staggered grid formulace

To je pro výuku v pořádku, ale pro robustnější CFD by se běžně použila přesnější a stabilnější formulace.

## Síť a diskretizace

Oblast:

- `x in [0, 1]`
- `y in [0, 1]`

Síť:

- `nx x ny` uzlů
- uniformní rozteč
- `dx = 1 / (nx - 1)`
- `dy = 1 / (ny - 1)`

Pole jsou ukládána do `std::vector<double>` a 2D indexace je mapována do 1D přes funkci:

```cpp
int idx(int i, int j, int nx) {
    return j * nx + i;
}
```

To znamená:

- `i` je index ve směru `x`
- `j` je index ve směru `y`
- lineární index je `j * nx + i`

## Okrajové podmínky

### Rychlost

Horní víko:

- `u = U_lid`
- `v = 0`

Ostatní stěny:

- `u = 0`
- `v = 0`

To odpovídá no-slip podmínce na pevných stěnách a pohybujícímu se hornímu víku.

### Tlak

Použitá jednoduchá volba:

- `dp/dn = 0` na levé, pravé a spodní stěně
- `p = 0` na horní stěně

Komentář:

- tlak v nestlačitelném proudění je určen jen do aditivní konstanty
- proto je potřeba jednu referenční hodnotu fixovat
- zde je zvoleno `p = 0` na horní stěně jako jednoduchá a čitelná volba

Je to praktické edukativní řešení, ne nutně jediná možná volba.

## Uživatelské parametry

Všechny hlavní parametry jsou na začátku `main.cpp`:

```cpp
// USER PARAMETERS
// Change these values first when you want to modify the case.
constexpr int nx = 41;
constexpr int ny = 41;
constexpr double dt = 0.001;
constexpr int nt = 500;
constexpr int nit = 50;
constexpr double rho = 1.0;
constexpr double nu = 0.1;
constexpr double u_lid = 1.0;
```

Význam:

- `nx`, `ny` - počet uzlů sítě
- `dt` - velikost časového kroku
- `nt` - počet časových kroků
- `nit` - počet iterací Poissonovy rovnice pro tlak v každém časovém kroku
- `rho` - hustota
- `nu` - kinematická viskozita
- `u_lid` - rychlost horního víka

## Jak projekt spustit

Přejdi do složky:

```bash
cd /home/ivand/projects/learning_cpp/cfd/programming_cfd_cases/01_cavity_case
```

Spusť překlad a výpočet:

```bash
./run.sh
```

Skript dělá:

```bash
g++ -std=c++20 -O2 -Wall -Wextra -pedantic main.cpp -o cavity
./cavity
```

Pokud chceš překlad spustit ručně:

```bash
g++ -std=c++20 -O2 -Wall -Wextra -pedantic main.cpp -o cavity
./cavity
```

## Jak změnit síť a fyzikální veličiny

Otevři `main.cpp` a uprav blok `USER PARAMETERS`.

### Příklad 1: jemnější síť

```cpp
constexpr int nx = 81;
constexpr int ny = 81;
```

### Příklad 2: delší simulace

```cpp
constexpr int nt = 2000;
```

### Příklad 3: menší viskozita

```cpp
constexpr double nu = 0.01;
```

### Příklad 4: rychlejší horní víko

```cpp
constexpr double u_lid = 2.0;
```

Po změně znovu spusť:

```bash
./run.sh
```

Pokud chceš změnit Reynoldsovo číslo:

- zvyš `u_lid`
- nebo sniž `nu`

Například:

```cpp
constexpr double nu = 0.01;
```

pak přibližně dostaneš:

```text
Re = 100
```

## Co čekat při změně parametrů

### Když zvýšíš `nx`, `ny`

- síť bude jemnější
- výsledek bude detailnější
- výpočet bude pomalejší
- tlakový solver bude potřebovat více práce

### Když zmenšíš `nu`

- proudění bude méně tlumené
- recirkulace bude výraznější
- explicitní schéma může být citlivější na stabilitu
- často je vhodné zmenšit `dt`

### Když zvětšíš `u_lid`

- proudění bude rychlejší
- konvekce bude silnější
- může být potřeba menší `dt`

### Když zvětšíš `dt`

- výpočet doběhne rychleji v menším počtu kroků
- ale může se zhoršit stabilita i přesnost
- při příliš velkém `dt` může řešení začít oscilovat nebo divergovat

### Když zvětšíš `nit`

- tlaková Poissonova rovnice se v každém kroku vyřeší lépe
- divergence bývá nižší
- výpočet bude pomalejší

## Konzolový výstup

Program během běhu vypisuje například:

```text
step 101/500  max|u|=1  max|v|=0.175738  div_l2=0.382712
```

Význam:

- `step` - aktuální časový krok
- `max|u|` - maximální absolutní hodnota složky rychlosti `u`
- `max|v|` - maximální absolutní hodnota složky rychlosti `v`
- `div_l2` - orientační norma divergence

Poznámka:

- divergence zde slouží hlavně jako jednoduchý indikátor kvality projekce
- není to plně rozvinuté reziduální monitorování

## Výstupní soubory

### CSV

Po doběhu vzniknou:

- `u.csv`
- `v.csv`
- `p.csv`

Každý soubor obsahuje hodnoty na mřížce ve formátu:

- řádky odpovídají směru `y`
- sloupce odpovídají směru `x`

### VTK pro ParaView

Soubor:

- `cavity.vtk`

Obsahuje:

- scalar pole `pressure`
- vector pole `velocity`

Použit byl jednoduchý formát `VTK legacy ASCII`, který ParaView otevře přímo.

## Jak otevřít výsledek v ParaView

Spuštění:

```bash
cd /home/ivand/projects/learning_cpp/cfd/programming_cfd_cases/01_cavity_case
paraview cavity.vtk
```

Nebo:

1. otevři `ParaView`
2. `File -> Open`
3. vyber `cavity.vtk`
4. klikni na `Apply`

## Doporučené zobrazení v ParaView

### Tlakové pole

- vyber pole `pressure`
- zobraz ho jako barevnou mapu

### Vektory rychlosti

- použij filtr `Glyph`
- jako vstupní vektor zvol `velocity`
- uprav hustotu a velikost šipek

### Proudnicové čáry

- použij filtr `Stream Tracer`
- jako vektorové pole zvol `velocity`

### Velikost rychlosti

Pokud chceš, můžeš v ParaView dopočítat velikost vektoru:

- použij `Calculator`
- výraz může být například `mag(velocity)`

## Typický workflow

Jednoduchý pracovní postup může být:

1. upravit parametry na začátku `main.cpp`
2. spustit `./run.sh`
3. otevřít `cavity.vtk` v ParaView
4. zkontrolovat tlak, vektory rychlosti a proudnice
5. podle potřeby upravit `dt`, `nit`, `nx`, `ny`, `nu`

## Důležité poznámky k numerice

Tento kód je záměrně jednoduchý. To znamená, že:

- není optimalizovaný na výkon
- není navržený pro vysoká Reynoldsova čísla
- nepoužívá pokročilé schéma pro konvekci
- nepoužívá staggered grid
- nemá sofistikované reziduální ani CFL řízení

Je vhodný hlavně pro:

- studium základní struktury CFD solveru
- první experimenty se sítí a parametry
- pochopení role tlaku v nestlačitelném proudění

## Troubleshooting

### Program se sice spustí, ale řešení vypadá nestabilně

Zkus:

- zmenšit `dt`
- zvýšit `nit`
- vrátit se na menší Reynoldsovo číslo

Typické první pokusy:

```cpp
constexpr double dt = 0.0005;
constexpr int nit = 100;
```

### `div_l2` je příliš vysoké

Tohle je u jednoduchého edukativního solveru do určité míry očekávatelné, ale může pomoct:

- zvýšit `nit`
- zjemnit síť
- zmenšit `dt`

### Výpočet je příliš pomalý

Zkus:

- menší síť, například `41 x 41`
- menší `nit`
- kratší simulaci přes nižší `nt`

### V ParaView nic nevidím

Zkontroluj:

- že jsi kliknul na `Apply`
- že máš otevřený `cavity.vtk`
- že je aktivní pole `pressure` nebo `velocity`
- že používáš vhodný typ zobrazení nebo filtr

### Šipky rychlosti jsou moc husté

V `Glyph` filtru:

- sniž počet vzorkovaných bodů
- zmenši scale factor
- případně nejdřív použij `Mask Points`

### Chci ostřejší víry, ale solver začne být citlivý

To bývá přirozené při:

- menší viskozitě
- větším `u_lid`
- jemnější dynamice proudění

První pomoc:

- sniž `dt`
- zvyš `nit`
- přecházej na náročnější parametry postupně

## Doporučené první experimenty

1. Změň síť z `41 x 41` na `81 x 81`.
2. Porovnej výsledné `cavity.vtk`.
3. Zkus změnit `nu` z `0.1` na `0.01`.
4. Sleduj, jak se mění recirkulace a stabilita.
5. Zkus zvýšit `nit` a pozoruj vliv na `div_l2`.

## Možná další rozšíření

### 1. Residual monitoring

Přidat:

- residuum tlakové Poissonovy rovnice
- residuum rychlosti mezi kroky
- jednodušší kritérium konvergence

### 2. Lepší export

Přidat:

- velikost rychlosti do VTK
- export v `VTU`
- více časových snapshotů

### 3. Numerická vylepšení

Přidat:

- staggered grid verzi
- lepší diskretizaci konvekce
- adaptivní `dt` podle CFL
- robustnější tlakový solver

## Rychlé minimum

Pokud si chceš zapamatovat jen to hlavní:

1. Parametry měníš na začátku `main.cpp`.
2. Spouštíš přes `./run.sh`.
3. Výsledek do ParaView je `cavity.vtk`.
