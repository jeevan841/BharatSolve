NAME          CRUDEBLEND
* Toy refinery blending problem for the SIH26119 demo.
* Blend LightCrude, HeavyCrude, Naphtha and Reformate into a gasoline pool
* meeting a volume target, a sulfur ceiling and an octane floor, at minimum
* cost. Weighted-average quality specs are written as linear constraints in
* the usual refinery-LP way: e.g. sulfur <= 0.5% of total volume becomes
* sum (sulfur_i - 0.5) * x_i <= 0.
ROWS
 N  COST
 G  VOLUME
 L  SULFUR
 G  OCTANE
 L  AVAIL1
 L  AVAIL2
 L  AVAIL3
 L  AVAIL4
COLUMNS
    LIGHTCR   COST            62.0   VOLUME           1.0
    LIGHTCR   SULFUR          -0.2   OCTANE            3.0
    LIGHTCR   AVAIL1           1.0
    HEAVYCR   COST            54.0   VOLUME           1.0
    HEAVYCR   SULFUR           1.3   OCTANE           -5.0
    HEAVYCR   AVAIL2           1.0
    NAPHTHA   COST            70.0   VOLUME           1.0
    NAPHTHA   SULFUR         -0.45   OCTANE            8.0
    NAPHTHA   AVAIL3           1.0
    REFORM    COST            58.0   VOLUME           1.0
    REFORM    SULFUR         -0.48   OCTANE           13.0
    REFORM    AVAIL4           1.0
RHS
    RHS       VOLUME       10000.0   SULFUR              0.0
    RHS       OCTANE           0.0   AVAIL1           6000.0
    RHS       AVAIL2        5000.0   AVAIL3           3000.0
    RHS       AVAIL4        4000.0
BOUNDS
ENDATA
