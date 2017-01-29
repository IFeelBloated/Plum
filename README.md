# Plum
©2017 IFeelBloated, Plum Python Module for VapourSynth

## License
LGPL v3.0

## Description
Plum is a sharpening/blind deconvolution suite with certain advanced features like Non-Local error, Block Matching, etc..

## Requirements
- [NNEDI3](https://github.com/dubhater/vapoursynth-nnedi3)
- [KNLMeansCL](https://github.com/Khanattila/KNLMeansCL)
- [BM3D](https://github.com/HomeOfVapourSynthEvolution/VapourSynth-BM3D)
- [FMTConv](https://github.com/EleonoreMizo/fmtconv)
- [VCFreq](http://www.avisynth.nl/users/vcmohan/vcfreq/vcfreq.html)
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
- RGB input will be converted to an opponent color space(YUV alike) and only luma will be processed still
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
The basic estimation performs sharpening with spatial self similarity to cancel out ringings, the sharpening kernel is not unsharpen masking, it's a blind deconvolution filter (assuming PSF is a circle).

workflow:
- do the blind deconvolution to enhance the overall sharpness.
- filter the result of deconvolution with Non-Local Errors, which shifts high frequency components from the deconvolved clip to the source clip without any new ringing (could enhance the existing ringing).
- clamp the result with a Local Errors unsharpen masking since most videos are not completely ringing-free, this makes sure the existing ringing won't be enhanced (at least won't be enhanced much).
- shrink the result down by 1 pixel
- repeat all the steps above a few times
- do another Non-Local Errors filtering to remove all residual ringing. (ringing inherited from the source clip)
- apply a cutoff filter to restore low frequency components.

```python
Basic(src, strength=6.4, a=32, h=64.0, radius=1, wn=0.48, scale=0.28, cutoff=32)
```
- strength<br />
  controls the iterating process, repeat floor(strength) and ceil(strength) times and blend them according to the fractional part of the strength.
- a<br />
  window size of the non-local error filtering.
- h<br />
  strength of the local errors filtering, greater value = more relaxed clamping.
- radius<br />
  radius of the deconvolution filter
- wn, scale<br />
  refer to VCFreq doc for more details.
- cutoff<br />
  strength of the cutoff filter, ranges from 0 (no low frequency protection) to 100 (almost no filtering)
