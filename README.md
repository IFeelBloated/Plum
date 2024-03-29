﻿# Plum
©2021 IFeelBloated, Plum Python Module for VapourSynth

## License
LGPL v3.0

## Description
Plum is a sharpening/blind deconvolution suite with certain advanced features like Non-Local error, Block Matching, etc..

## Requirements
- [NNEDI3](https://github.com/dubhater/vapoursynth-nnedi3)
- [KNLMeansCL](https://github.com/Khanattila/KNLMeansCL)
- [BM3D](https://github.com/HomeOfVapourSynthEvolution/VapourSynth-BM3D)
- [FMTConv](https://github.com/EleonoreMizo/fmtconv)
- [VCM](http://www.avisynth.nl/users/vcmohan/vcm/vcm.html)
- [MVTools (floating point ver)](https://github.com/IFeelBloated/vapoursynth-mvtools-sf/tree/master)

## Function List
- Super
- Basic
- Final

## Formats
- Bit Depth: 32bits floating point
- Color Space: Gray, RGB, YUV4XXPS
- Scan Type: Progressive

## Notes
- Only Y will be processed when the input is YUV, UV will be simply copied from the input
- RGB input will be converted to an opponent color space(YUV alike) and only luma will be processed
- **NO** scene change policy provided, take [Wobbly](https://github.com/dubhater/Wobbly) and cut each scene out and process them individually
- **QUALITY**: cutting edge
- **PERFORMANCE**: close to abysmal

## Details
### Super
Optional, it helps improve the precision of sub-pixel motion estimation and compensation, use it and get a quality boost or don't and get a performance boost
```python
Super(src, pel=4)
```
- src<br />
  clip to be processed
- pel<br />
  sub-pixel precision, could be 2 or 4, 2 = precision by half a pixel, 4 = precision by quarter a pixel.

### Basic
The basic estimation performs sharpening with spatial self similarity to cancel out ringing, the sharpening kernel is not unsharpen masking, it's a blind deconvolution filter (assuming PSF is a circle).

workflow:
- do the blind deconvolution to enhance the overall sharpness.
- filter the result of deconvolution with supersampled Non-Local Errors, which shifts high frequency components from the deconvolved clip to the source clip without any new ringing/aliasing (could enhance the existing ringing/aliasing).
- clamp the result with a Local Errors unsharpen masking since most videos are not completely ringing-free, this makes sure the existing ringing won't be enhanced (at least won't be enhanced much).
- shrink the result down by 1 pixel
- repeat all steps above a few times
- do another Non-Local Errors filtering to remove all residual ringing. (ringing inherited from the source clip)
- apply a cutoff filter to restore low frequency components.

```python
Basic(src, strength=3.20, a=32, h=[6.4, 64.0], radius=1, wn=0.48, scale=0.28, cutoff=24)
```
- strength<br />
  controls the iterating process, repeat floor(strength) and ceil(strength) times and blend them according to the fractional part of the strength.
- a<br />
  window size of the non-local error filtering.
- h<br />
  h[0]: standard deviation of the non-local error filtering, default value is pretty balanced.<br />
  h[1]: strength of the local errors filtering, greater value = more relaxed clamping.
- radius<br />
  radius of the deconvolution filter
- wn, scale<br />
  refer to VCFreq doc for more details.
- cutoff<br />
  strength of the cutoff filter, ranges from 0 (no low frequency protection) to 100 (almost no filtering)

### Final
The final estimation adjusts the basic estimation using temporal self similarity to cancel out noise and residual aliasing.

workflow:
- do a motion compensated temporal averaging to the difference between the basic estimation and the source clip, then apply the stabilized difference back to the source clip.
- clamp the result temporally with motion compensation.
- get the difference between the result and the source clip and amplify it non-linearly, then apply it back.
- apply a cutoff filter to restore low frequency components.

```python
Final(src, super=[None, None], radius=6, pel=4, sad=400.0, flexibility=0.64, strength=3.20, constants=[1.49, 1.272, None], cutoff=12, freq_margin=20)
```
- super<br />
  optional, clips generated by Plum.Super
- radius<br />
  temporal radius, frames that fall in [current frame - radius, current frame + radius] will be referenced
- sad<br />
  SAD threshold of the motion compensation, refer to MVTools doc for more details
- flexibility<br />
  flexibility of the motion compensated temporal clamping, on a scale of [0.0, 1.0], greater value = more relaxed clamping
- strength<br />
  general amplitude of the non-linear amplification function.
- constants<br />
  parameters related to the non-linear amplification function<br />
  constants[0]: modifier for the amplification function<br />
  constants[1]: exponent for the amplification function<br />
  constants[2]: suppression to the very small differences, default = strength + 0.1
- freq_margin<br />
  frequency margin of the cutoff threshold, larger value = more delicate and less agressive sharpening

## Demos
- A
```python
ref = Plum.Basic(clip, strength=6.4, cutoff=32)
clip = Plum.Final([clip, ref], [Plum.Super(clip), Plum.Super(ref)], strength=1.8, freq_margin=12)
```
![](http://i.imgur.com/X3mG5NX.png)
![](http://i.imgur.com/CMwOlYx.png)
- B
```python
ref = Plum.Basic(clip)
clip = Plum.Final([clip, ref], [Plum.Super(clip), Plum.Super(ref)], cutoff=8, freq_margin=12)
```
![](http://i.imgur.com/8PqTPbC.png)
![](http://i.imgur.com/E5Zi0TO.png)
- C
```python
ref = Plum.Basic(clip)
clip = Plum.Final([clip, ref], [Plum.Super(clip), Plum.Super(ref)])
```
![](http://i.imgur.com/JdbJPuM.png)
![](http://i.imgur.com/BwfCu6E.png)
