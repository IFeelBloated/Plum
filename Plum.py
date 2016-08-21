import vapoursynth as vs
import mvmulti

fmtc_args                 = dict(fulls=True, fulld=True)
msuper_args               = dict(hpad=32, vpad=32, sharp=2, levels=0, chroma=False)
manalyze_args             = dict(search=3, truemotion=False, trymany=True, levels=0, badrange=-24, divide=0, dct=0, chroma=False)
mrecalculate_args         = dict(truemotion=False, search=3, smooth=1, divide=0, dct=0, chroma=False)
mdegrain_args             = dict(plane=0)
conv_args                 = dict(matrix=[1, 2, 1, 2, 4, 2, 1, 2, 1])
deconv_args               = dict(line=0, wn=0.48, fr=25, scale=0.28)
dfttest_args              = dict(smode=0, sosize=0, tbsize=1, tosize=0, tmode=0)
nnedi_args                = dict(field=1, dh=True, nns=4, qual=2, etype=1, nsize=0)

class helpers:
      def cutoff(low, hi, p):
          core            = vs.get_core()
          Resample        = core.fmtc.resample
          MakeDiff        = core.std.MakeDiff
          MergeDiff       = core.std.MergeDiff
          def inline(src):
              upsmp       = Resample(src, src.width * 2, src.height * 2, kernel="gauss", a1=100, **fmtc_args)
              clip        = Resample(upsmp, src.width, src.height, kernel="gauss", a1=p, **fmtc_args)
              return clip
          hif             = MakeDiff(hi, inline(hi))
          clip            = MergeDiff(inline(low), hif)
          return clip
      def padding(src, left, right, top, bottom):
          core            = vs.get_core()
          Resample        = core.fmtc.resample
          w               = src.width
          h               = src.height
          clip            = Resample(src, w+left+right, h+top+bottom, -left, -top, w+left+right, h+top+bottom, kernel="point", **fmtc_args)
          return clip
      def deconvolution(src, radius):
          core            = vs.get_core()
          FQSharp         = core.vcfreq.Sharp
          sharp           = FQSharp(src, x=radius, y=radius, **deconv_args)
          clip            = helpers.cutoff(src, sharp, 1)
          return clip
      def convolution(src, strength):
          core            = vs.get_core()
          Resample        = core.fmtc.resample
          NNEDI           = core.nnedi3.nnedi3
          Transpose       = core.std.Transpose
          MakeDiff        = core.std.MakeDiff
          MergeDiff       = core.std.MergeDiff
          w               = src.width
          h               = src.height
          supersampled    = Transpose(NNEDI(Transpose(NNEDI(src, **nnedi_args)), **nnedi_args))
          blur            = Resample(supersampled, w*8, h*8, kernel="cubic", a1=strength, a2=0, **fmtc_args)
          sharp           = Resample(supersampled, w*8, h*8, kernel="cubic", a1=-strength, a2=0, **fmtc_args)
          dif             = Resample(MakeDiff(sharp, blur), w, h, sx=-0.5, sy=-0.5, kernel="cubic", a1=-1, a2=0, **fmtc_args)
          clip            = MergeDiff(src, dif)
          return clip
      def shrink(src):
          core            = vs.get_core()
          Convolution     = core.std.Convolution
          Expr            = core.std.Expr
          Crop            = core.std.CropRel
          MakeDiff        = core.std.MakeDiff
          Median          = core.std.Median
          MergeDiff       = core.std.MergeDiff
          blur            = Median(src)
          dif             = MakeDiff(blur, src)
          convD           = Convolution(dif, **conv_args)
          DD              = MakeDiff(dif, convD)
          convDD          = Convolution(DD, **conv_args)
          DDD             = Expr([DD, convDD], ["x y - x 0.5 - * 0 < 0.5 x y - abs x 0.5 - abs < x y - 0.5 + x ? ?"])
          dif             = MakeDiff(dif, DDD)
          convD           = Convolution(dif, **conv_args)
          dif             = Expr([dif, convD], ["y 0.5 - abs x 0.5 - abs > y 0.5 ?"])
          clip            = MergeDiff(src, dif)
          return clip
      def nlerror(src, a, h, ref):
          core            = vs.get_core()
          Crop            = core.std.CropRel
          KNLMeansCL      = core.knlm.KNLMeansCL
          pad             = helpers.padding(src, a, a, a, a)
          ref             = helpers.padding(ref, a, a, a, a)
          nlm             = KNLMeansCL(pad, d=0, a=a, s=0, h=h, rclip=ref)
          clip            = Crop(nlm, a, a, a, a)
          return clip
      def extremum_multi(src, radius, mode):
          core            = vs.get_core()
          SelectEvery     = core.std.SelectEvery
          Expr            = core.std.Expr
          clip            = SelectEvery(src, radius * 2 + 1, 0)
          for i in range(1, radius * 2 + 1):
              clip        = Expr([clip, SelectEvery(src, radius * 2 + 1, i)], "x y " + mode)
          return clip
      def clamp(src, bright_limit, dark_limit, overshoot, undershoot):
          core            = vs.get_core()
          Expr            = core.std.Expr
          clip            = Expr([src, bright_limit, dark_limit], ["x y {os} + > y {os} + x ? z {us} - < z {us} - x ?".format(os=overshoot, us=undershoot)])
          return clip

class internal:
      def super(src, pel):
          core            = vs.get_core()
          NNEDI           = core.nnedi3.nnedi3
          MakeDiff        = core.std.MakeDiff
          MergeDiff       = core.std.MergeDiff
          Transpose       = core.std.Transpose
          clip            = [None, None]
          clip[0]         = Transpose(NNEDI(Transpose(NNEDI(src[0], **nnedi_args)), **nnedi_args))
          if pel == 4:
             clip[0]      = Transpose(NNEDI(Transpose(NNEDI(clip[0], **nnedi_args)), **nnedi_args))
          if src[1] is not None:
             src[1]       = MergeDiff(src[0], src[1])
             clip[1]      = Transpose(NNEDI(Transpose(NNEDI(src[1], **nnedi_args)), **nnedi_args))
             if pel == 4:
                clip[1]   = Transpose(NNEDI(Transpose(NNEDI(clip[1], **nnedi_args)), **nnedi_args))
             clip[0]      = MakeDiff(clip[1], clip[0])
          return clip[0]
      def basic(src, iterate, a, h, deconv_radius, conv_strength, mode):
          core            = vs.get_core()
          Expr            = core.std.Expr
          MakeDiff        = core.std.MakeDiff
          MergeDiff       = core.std.MergeDiff
          if mode == "deconvolution":
             sharp        = helpers.deconvolution(src, deconv_radius)
          else:
             sharp        = helpers.convolution(src, conv_strength)
          sharp           = helpers.nlerror(src, a[0], 0.001, sharp)
          local_error     = helpers.nlerror(src, a[1], h, src)
          local_limit     = MergeDiff(src, MakeDiff(src, local_error))
          limited         = Expr([sharp, local_limit, src], ["x z - abs y z - abs > y x ?"])
          clip            = helpers.shrink(limited)
          iterate        -= 1
          if iterate == 0:
             return clip
          else:
             return internal.basic(clip, iterate, a, h, deconv_radius, conv_strength, mode)
      def final(src, super, radius, pel, sad, strength, constants, attenuate_window, attenuate, cutoff):
          core            = vs.get_core()
          DFTTest         = core.dfttest.DFTTest
          MSuper          = core.mvsf.Super
          MAnalyze        = mvmulti.Analyze
          MRecalculate    = mvmulti.Recalculate
          MDegrainN       = mvmulti.DegrainN
          MCompensate     = mvmulti.Compensate
          Expr            = core.std.Expr
          MakeDiff        = core.std.MakeDiff
          MergeDiff       = core.std.MergeDiff
          expression      = "{x} {y} - abs {lstr} / 1 {pstr} / pow {sstr} * {x} {y} - {x} {y} - abs 0.001 + / * {x} {y} - 2 pow {x} {y} - 2 pow {ldmp} + / * 256 / y +".format(lstr=constants[0], pstr=constants[1], sstr=strength, ldmp=constants[2], x="x 256 *", y="y 256 *")
          blankdif        = Expr(src[0], "0.5")
          supersoft       = MSuper(src[0], pelclip=super[0], rfilter=4, pel=pel, **msuper_args)
          supersharp      = MSuper(src[0], pelclip=super[0], rfilter=2, pel=pel, **msuper_args)
          superdif        = MSuper(src[1], pelclip=super[1], rfilter=2, pel=pel, **msuper_args)
          superlimit      = MSuper(src[2], pelclip=super[2], rfilter=2, pel=pel, **msuper_args)
          vmulti          = MAnalyze(supersoft, tr=radius, overlap=16, blksize=32, **manalyze_args)
          vmulti          = MRecalculate(supersoft, vmulti, tr=radius, overlap=8, blksize=16, thsad=sad/2.0, **mrecalculate_args)
          vmulti          = MRecalculate(supersharp, vmulti, tr=radius, overlap=4, blksize=8, thsad=sad/2.0, **mrecalculate_args)
          vmulti          = MRecalculate(supersharp, vmulti, tr=radius, overlap=2, blksize=4, thsad=sad/2.0, **mrecalculate_args)
          averaged_dif    = MDegrainN(blankdif, superdif, vmulti, tr=radius, thsad=10000.0, thscd1=10000.0, thscd2=255.0, **mdegrain_args)
          compensated     = MCompensate(src[0], superlimit, vmulti, tr=radius, thsad=sad, thscd1=10000.0, thscd2=255.0)
          bright_limit    = helpers.extremum_multi(compensated, radius, "max")
          dark_limit      = helpers.extremum_multi(compensated, radius, "min")
          averaged        = MergeDiff(src[0], averaged_dif)
          clamped         = helpers.clamp(averaged, bright_limit, dark_limit, 0.0, 0.0)
          amplified       = Expr([clamped, src[0]], expression)
          dif             = MakeDiff(amplified, src[0])
          dif             = DFTTest(dif, sbsize=attenuate_window, sstring=attenuate, **dfttest_args)
          sharp           = MergeDiff(src[0], dif)
          clip            = helpers.cutoff(src[0], sharp, cutoff)
          return clip

def Super(src, pel=4):
    core                  = vs.get_core()
    RGB2OPP               = core.bm3d.RGB2OPP
    ShufflePlanes         = core.std.ShufflePlanes
    SetFieldBased         = core.std.SetFieldBased
    if not isinstance(src, list):
       raise TypeError("Plum.Super: src has to be an array!")
    elif len(src) != 2:
       raise RuntimeError("Plum.Super: src has to contain 2 elements exactly!")
    if not isinstance(src[0], vs.VideoNode):
       raise TypeError("Plum.Super: src[0] has to be a video clip!")
    elif src[0].format.sample_type != vs.FLOAT or src[0].format.bits_per_sample < 32:
       raise TypeError("Plum.Super: the sample type of src[0] has to be single precision!")
    if not isinstance(src[1], vs.VideoNode) and src[1] is not None:
       raise TypeError("Plum.Super: src[1] has to be a video clip or None!")
    elif src[1] is not None:
       if src[1].format.id != vs.GRAYS:
          raise RuntimeError("Plum.Super: corrupted basic estimation!")
    if not isinstance(pel, int):
       raise TypeError("Plum.Super: pel has to be an integer!")
    elif pel != 2 and pel != 4:
       raise RuntimeError("Plum.Super: pel has to be 2 or 4!")
    src[0]                = SetFieldBased(src[0], 0)
    src[1]                = SetFieldBased(src[1], 0) if src[1] is not None else None
    colorspace            = src[0].format.color_family
    if colorspace == vs.RGB:
       src[0]             = RGB2OPP(src[0], 1)
    if colorspace != vs.GRAY:
       src[0]             = ShufflePlanes(src[0], 0, vs.GRAY)
    clip                  = internal.super(src, pel)
    return clip

def Basic(src, iterate=3, a=[32, 1], h=64.0, deconv_radius=1, conv_strength=3.2, mode="deconvolution"):
    core                  = vs.get_core()
    RGB2OPP               = core.bm3d.RGB2OPP
    MakeDiff              = core.std.MakeDiff
    ShufflePlanes         = core.std.ShufflePlanes
    SetFieldBased         = core.std.SetFieldBased
    if not isinstance(src, vs.VideoNode):
       raise TypeError("Plum.Basic: src has to be a video clip!")
    elif src.format.sample_type != vs.FLOAT or src.format.bits_per_sample < 32:
       raise TypeError("Plum.Basic: the sample type of src has to be single precision!")
    if not isinstance(iterate, int):
       raise TypeError("Plum.Basic: iterate has to be an integer!")
    elif iterate < 1:
       raise RuntimeError("Plum.Basic: iterate has to be greater than 0!")
    if not isinstance(a, list):
       raise TypeError("Plum.Basic: a has to be an array!")
    elif len(a) != 2:
       raise RuntimeError("Plum.Basic: a has to contain 2 elements exactly!")
    elif not isinstance(a[0], int) or not isinstance(a[1], int):
       raise TypeError("Plum.Basic: elements in a must be integers!")
    if not isinstance(h, float) and not isinstance(h, int):
       raise TypeError("Plum.Basic: h has to be a real number!")
    elif h <= 0:
       raise RuntimeError("Plum.Basic: h has to be greater than 0!")
    if not isinstance(deconv_radius, int):
       raise TypeError("Plum.Basic: deconv_radius has to be an integer!")
    elif deconv_radius < 1:
       raise RuntimeError("Plum.Basic: deconv_radius has to be greater than 0!")
    if not isinstance(conv_strength, float) and not isinstance(conv_strength, int):
       raise TypeError("Plum.Basic: conv_strength has to be a real number!")
    elif conv_strength <= 0:
       raise RuntimeError("Plum.Basic: conv_strength has to be greater than 0!")
    if not isinstance(mode, str):
       raise TypeError("Plum.Basic: mode has to be a string!")
    elif mode.lower() != "deconvolution" and mode.lower() != "convolution":
       raise NotImplementedError("Plum.Basic: Undefined mode!")
    src                   = SetFieldBased(src, 0)
    colorspace            = src.format.color_family
    if colorspace == vs.RGB:
       src                = RGB2OPP(src, 1)
    if colorspace != vs.GRAY:
       src                = ShufflePlanes(src, 0, vs.GRAY)
    clip                  = internal.basic(src, iterate, a, h, deconv_radius, conv_strength, mode.lower())
    if mode.lower() == "deconvolution":
       clip               = MakeDiff(clip, src)
    return clip

def Final(src, super=[None, None, None], radius=6, pel=4, sad=400.0, strength=1.64, constants=[1.49, 1.272, None], \
          attenuate_window=9, attenuate="0.0:1024.0 0.16:1024.0 0.32:1024.0 0.48:128.0 0.56:0.0 1.0:0.0", cutoff=12):
    core                  = vs.get_core()
    RGB2OPP               = core.bm3d.RGB2OPP
    OPP2RGB               = core.bm3d.OPP2RGB
    ShufflePlanes         = core.std.ShufflePlanes
    SetFieldBased         = core.std.SetFieldBased
    if not isinstance(src, list):
       raise TypeError("Plum.Final: src has to be an array!")
    elif len(src) != 3:
       raise RuntimeError("Plum.Final: src has to contain 3 elements exactly!")
    elif not isinstance(src[0], vs.VideoNode) or not isinstance(src[1], vs.VideoNode) or not isinstance(src[2], vs.VideoNode):
       raise TypeError("Plum.Final: elements in src must be video clips!")
    elif src[0].format.sample_type != vs.FLOAT or src[0].format.bits_per_sample < 32:
       raise TypeError("Plum.Final: the sample type of src[0] has to be single precision!")
    elif src[1].format.id != vs.GRAYS or src[2].format.id != vs.GRAYS:
       raise RuntimeError("Plum.Final: corrupted basic estimation!")
    if not isinstance(super, list):
       raise TypeError("Plum.Final: super has to be an array!")
    elif len(super) != 3:
       raise RuntimeError("Plum.Final: super has to contain 3 elements exactly!")
    for i in range(3):
        if not isinstance(super[i], vs.VideoNode) and super[i] is not None:
           raise TypeError("Plum.Final: elements in super must be video clips or None!")
        elif super[i] is not None:
           if super[i].format.id != vs.GRAYS:
              raise RuntimeError("Plum.Final: corrupted super clips!")
    if not isinstance(radius, int):
       raise TypeError("Plum.Final: radius has to be an integer!")
    elif radius < 1:
       raise RuntimeError("Plum.Final: radius has to be greater than 0!")
    if not isinstance(pel, int):
       raise TypeError("Plum.Final: pel has to be an integer!")
    elif pel != 1 and pel != 2 and pel != 4:
       raise RuntimeError("Plum.Final: pel has to be 1, 2 or 4!")
    if not isinstance(sad, float) and not isinstance(sad, int):
       raise TypeError("Plum.Final: sad has to be a real number!")
    elif sad <= 0:
       raise RuntimeError("Plum.Final: sad has to be greater than 0!")
    if not isinstance(strength, float) and not isinstance(strength, int):
       raise TypeError("Plum.Final: strength has to be a real number!")
    elif strength <= 0:
       raise RuntimeError("Plum.Final: strength has to be greater than 0!")
    if not isinstance(constants, list):
       raise TypeError("Plum.Final: constants parameter has to be an array!")
    elif len(constants) != 3:
       raise RuntimeError("Plum.Final: constants parameter has to contain 3 elements exactly!")
    for i in range(2):
        if not isinstance(constants[i], float) and not isinstance(constants[i], int):
           raise TypeError("Plum.Final: elements in constants must be real numbers!")
    if not isinstance(constants[2], float) and not isinstance(constants[2], int) and constants[2] is not None:
       raise TypeError("Plum.Final: constants[2] has to be a real number or None!")
    if not isinstance(attenuate_window, int):
       raise TypeError("Plum.Final: attenuate_window has to be an integer!")
    elif attenuate_window < 1 or attenuate_window % 2 != 1:
       raise RuntimeError("Plum.Final: attenuate_window has to be a positive odd integer!")
    if not isinstance(attenuate, str):
       raise TypeError("Plum.Final: attenuate has to be a string!")
    if not isinstance(cutoff, int):
       raise TypeError("Plum.Final: cutoff has to be an integer!")
    elif cutoff < 1 or cutoff > 100:
       raise RuntimeError("Plum.Final: cutoff must fall in(0, 100]!")
    constants[2]          = strength + 0.1 if constants[2] is None else constants[2]
    for i in range(3):
        src[i]            = SetFieldBased(src[i], 0)
        super[i]          = SetFieldBased(super[i], 0) if super[i] is not None else None
    colorspace            = src[0].format.color_family
    if colorspace == vs.RGB:
       src[0]             = RGB2OPP(src[0], 1)
    if colorspace != vs.GRAY:
       src_color          = src[0]
       src[0]             = ShufflePlanes(src[0], 0, vs.GRAY)
    clip                  = internal.final(src, super, radius, pel, sad, strength, constants, attenuate_window, attenuate, cutoff)
    if colorspace != vs.GRAY:
       clip               = ShufflePlanes([clip, src_color], [0, 1, 2], vs.YUV)
    if colorspace == vs.RGB:
       clip               = OPP2RGB(clip, 1)
    return clip
