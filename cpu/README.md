# FIMD-CPU: Fast Isolated Marker Detection for CPUs

This directory contains the CPU implementation of the FIMD algorithm. The universal implementation is located in the `template.c` file. This template is processed by the `generate.py` Python script to produce the final C implementation. 

CMake compiles two targets:

* `fimd_cpu` - A shared library exposing a detection function in the header file `fimd_cpu.h`. The radii must be specified in the file `CMakelists.txt` before compilation.
* `fimd_cpu_example` - An executable for testing the detection function with source code in the `example.c` file.

## Circle boundary and interior generation (example)
The boundary and interior points are generated by the Python script in the final evaluation order. Below is an example of verbose output for a radius of 6:

```txt
Generated Bresenham circle:
-- Boundary points: 32
-- Interior points: 89
Visualization:
                    B10 B02 B06
            B26 B18 I55 I21 I53 B14 B22
        B32 I87 I75 I47 I17 I45 I73 I85 B30
    B28 I89 I81 I67 I39 I13 I37 I65 I79 I88 B24
    B20 I77 I69 I61 I31 I09 I29 I59 I68 I76 B16
B12 I57 I49 I41 I33 I25 I05 I23 I32 I40 I48 I56 B08
B03 I20 I16 I12 I08 I04 I01 I03 I07 I11 I15 I19 B04
B07 I52 I44 I36 I28 I24 I02 I22 I27 I35 I43 I51 B11
    B15 I72 I64 I60 I30 I06 I26 I58 I63 I71 B19
    B23 I84 I80 I66 I38 I10 I34 I62 I78 I83 B27
        B29 I86 I74 I46 I14 I42 I70 I82 B31
            B21 B13 I54 I18 I50 B17 B25
                    B05 B01 B09

Pixel evaluation order:
-- Boundary points: 32
-- Interior points: 45
Visualization:
                    B18 B02 B32
            B22 B10             B16 B28
        B06                             B08
    B26                                     B24
    B14                                     B12
B30                                             B20
B04                     I01 I02 I03 I04 I05 I06 B03
B19 I07 I08 I09 I10 I11 I12 I13 I14 I15 I16 I17 B29
    B11 I18 I19 I20 I21 I22 I23 I24 I25 I26 B13
    B23 I27 I28 I29 I30 I31 I32 I33 I34 I35 B25
        B07 I36 I37 I38 I39 I40 I41 I42 B05
            B27 B15 I43 I44 I45 B09 B21
                    B31 B01 B17
```

## Detection output (example)
Detection output with (x, y) coordinates for a dataset sample `1619240573769481609.bin` and selected radii 1 to 10:

```txt
FIMD-CPU r=1: detected 1 markers, 2591 sun points.
Markers: [(273,318)]

FIMD-CPU r=2: detected 12 markers, 364 sun points.
Markers: [(266,194),(257,318),(283,325),(221,354),(229,358),(520,382),(523,382),(515,384),(411,407),(379,415),(386,416),(374,416)]

FIMD-CPU r=3: detected 11 markers, 181 sun points.
Markers: [(266,194),(272,316),(257,317),(283,324),(221,354),(229,358),(520,382),(515,384),(379,415),(386,416),(374,416)]

FIMD-CPU r=4: detected 12 markers, 119 sun points.
Markers: [(266,194),(272,315),(256,317),(283,324),(221,354),(229,358),(520,382),(515,384),(412,407),(379,415),(386,416),(374,416)]

FIMD-CPU r=5: detected 11 markers, 82 sun points.
Markers: [(266,194),(272,315),(256,317),(283,324),(221,354),(229,358),(520,382),(515,384),(411,407),(379,415),(386,416)]

FIMD-CPU r=6: detected 10 markers, 52 sun points.
Markers: [(266,194),(272,315),(256,317),(283,324),(221,354),(229,358),(520,382),(411,407),(379,415),(386,416)]

FIMD-CPU r=7: detected 9 markers, 41 sun points.
Markers: [(266,194),(272,315),(256,317),(283,324),(221,354),(229,358),(520,382),(411,407),(379,415)]

FIMD-CPU r=8: detected 9 markers, 34 sun points.
Markers: [(266,194),(272,315),(256,317),(283,324),(221,354),(229,358),(520,382),(411,407),(379,415)]

FIMD-CPU r=9: detected 8 markers, 24 sun points.
Markers: [(266,194),(272,315),(256,317),(283,324),(221,354),(520,382),(411,407),(379,415)]

FIMD-CPU r=10: detected 8 markers, 20 sun points.
Markers: [(266,194),(272,315),(256,317),(283,324),(221,354),(520,382),(411,407),(379,415)]
```

## Symbol table for the shared library (example)

Output of the `objdump -TC libfimd_cpu.so` command:
```txt
libfimd_cpu.so:     file format elf64-x86-64

DYNAMIC SYMBOL TABLE:
0000000000000000      DF *UND*    0000000000000000  GLIBC_2.2.5 free
0000000000000000  w   D  *UND*    0000000000000000              _ITM_deregisterTMCloneTable
0000000000000000  w   D  *UND*    0000000000000000              __gmon_start__
0000000000000000      DF *UND*    0000000000000000  GLIBC_2.14  memcpy
0000000000000000      DF *UND*    0000000000000000  GLIBC_2.2.5 malloc
0000000000000000  w   D  *UND*    0000000000000000              _ITM_registerTMCloneTable
0000000000000000  w   DF *UND*    0000000000000000  GLIBC_2.2.5 __cxa_finalize
0000000000001a60 g    DF .text    0000000000000534  Base        fimd_r3
00000000000011d0 g    DF .text    0000000000000320  Base        fimd_cpu_detect
0000000000001520 g    DF .text    0000000000000008  Base        fimd_cpu_get_radii
0000000000001fa0 g    DF .text    0000000000000704  Base        fimd_r4
00000000000026b0 g    DF .text    00000000000008e8  Base        fimd_r5
0000000000002fa0 g    DF .text    0000000000000c7e  Base        fimd_r6
0000000000003c20 g    DF .text    0000000000000e99  Base        fimd_r7
0000000000001500 g    DF .text    0000000000000006  Base        fimd_cpu_image_height
0000000000004ac0 g    DF .text    000000000000115f  Base        fimd_r8
0000000000005c20 g    DF .text    00000000000015a4  Base        fimd_r9
0000000000001530 g    DF .text    0000000000000006  Base        fimd_cpu_get_max_markers_count
0000000000001560 g    DF .text    0000000000000006  Base        fimd_cpu_get_threshold_sun
0000000000001570 g    DF .text    0000000000000006  Base        fimd_cpu_get_threshold_diff
0000000000001540 g    DF .text    0000000000000006  Base        fimd_cpu_get_max_sun_points_count
00000000000071d0 g    DF .text    00000000000018d7  Base        fimd_r10
0000000000001580 g    DF .text    0000000000000006  Base        fimd_cpu_get_termination_sequence
0000000000001550 g    DF .text    0000000000000006  Base        fimd_cpu_get_threshold_marker
0000000000009030 g    DO .rodata  0000000000000028  Base        fimd_radii_list
00000000000014f0 g    DF .text    0000000000000006  Base        fimd_cpu_image_width
0000000000001590 g    DF .text    0000000000000156  Base        fimd_r1
00000000000016f0 g    DF .text    0000000000000366  Base        fimd_r2
0000000000001510 g    DF .text    0000000000000006  Base        fimd_cpu_get_radii_count
```


## ARMv6 Assembly implementation for MCUs

The file `template_ARMv6.s` contains an ARMv6 assembly implementation of the FIMD algorithm, which is designed for low-power MCUs, such as STM32. The assembly template is not compiled by the CMake build system and must be manually processed by the generator Python script `generate.py`. The generated output is ready for integration into the target STM32CubeIDE project. A header file `fimd_armv6.h` (referenced by the template file) should contain the requisite macros with image dimensions, selected thresholds, radius, detection limits and termination sequence (similar to `CMakeLists.txt`).
