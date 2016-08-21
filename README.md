# Plum
©2016 IFeelBloated, Plum Python Module for VapourSynth

## License
LGPL v3.0

## Description
Plum is a sharpening/blind deconvolution suite with certain advanced features like Non-Local error, Block Matching, etc..

## Requirements
- [NNEDI3](https://github.com/dubhater/vapoursynth-nnedi3)
- [KNLMeansCL](https://github.com/Khanattila/KNLMeansCL)
- [BM3D](https://github.com/HomeOfVapourSynthEvolution/VapourSynth-BM3D)
- [DFTTest](https://github.com/HomeOfVapourSynthEvolution/VapourSynth-DFTTest)
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
- RGB input will be converted to an opponent color space(YUV alike) and only luminance will be processed still
- **QUALITY**: cutting edge
- **PERFORMANCE**: close to abysmal
