# 01 — First run

BBPA is a function object: a block inside `functions{}` in
`system/controlDict`. The smallest working dictionary is:

```c
functions
{
    bbpa
    {
        type            BBPA;
        libs            ("libPhaseAveraging.so");

        executeControl  timeStep;
        executeInterval 1;
        writeControl    writeTime;

        fields          (U);
        period          1.0;
    }
}
```

Six of these entries are structural. `type` and `libs` load the code
(chapter 03). The three scheduling entries make BBPA see every solver
step and write when the solver writes (chapter 04). `fields` says what
to average (chapter 05). The only physics you must supply is `period`:
the cycle length of your flow, in seconds.

Everything else has a default: 100 bins per cycle, no discarded cycles,
no cycle limit, phase origin aligned to the run start, left-aligned
windows, every bin accumulated and written, one dump per completed
cycle once `output{}` is given (chapter 10 has the exact rules), and
all statistical products on.

## Run it

The reference case is `tutorials/bbpaVerify`: 320 cells, a sinusoidal
inlet with period 1 s, 10 bins, 3 simulated seconds, about 10 s of wall
clock. Its `controlDict` is the annotated version of the block above.

```bash
source /opt/openfoam14/etc/bashrc
cd src/phaseAveraging && wmake libso && cd -     # build the library once

cd tutorials/bbpaVerify
./Allclean && ./Allrun
```

While it runs, the log shows BBPA's bookkeeping: the configuration
echo at startup, then per write a line like
`Cycles accumulated: 2` and `Bins written: 10/10`.

## Check it

```bash
python3 check_bbpa.py
```

This re-derives the cross-cycle statistics from BBPA's own per-cycle
output with independent numpy code and compares. It ends with `PASS`.
Chapter 02 explains what the files it read actually are.
