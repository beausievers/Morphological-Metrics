#
# ==Introduction
#
# This module contains implementations of the metrics described by Larry
# Polansky in his article _Morphological_ _Metrics_ (_MM_ in this document).
# The abstract of _MM_ introduces the idea of measuring musical distance:
#
# "This article describes a number of methods for measuring morphological
# similarity in music. A variety of _metric_ _equations_ are described which 
# attempt to understand certain kinds of musical and perceptual distances 
# between pairs of shapes in any musical parameter. Different techniques for 
# applying these metrics, such as weightings, variants on each class of 
# metric, and multidimensional and multimetric combinations are also 
# described." _MM_ p. 289
#
# _MM_ defines morphologies, or morphs, as ordered sets. This module treats 
# morphs as 1-dimensional vectors of numerical values, and depends on the
# NArray library for its vector implementation.
#
# The original article and other supporting information can be found on 
# Larry's website: http://eamusic.dartmouth.edu/~larry/mutationsFAQ.html
#
# In addition to the techniques described in _MM_, this module contains tools 
# for generating morphs and sets of morphs which satisfy certain constraints
# regarding their positions in metric space. Though these tools were designed
# with music composition and analysis in mind, they may find application in 
# other areas, such as visual art, computer graphics, psychological experiment 
# design, and the natural sciences. 
#
# ==Simple Usage Examples
# 
# Metrics can be called in a number of ways. The simplest is to prepend 
# "dist_" to the name of the metric you want to call, like so:
#
#  > m = NArray[1,2,3,4]
#  > n = NArray[1,3,5,7]
#  > MM.dist_ocm(m,n)
#  => 0.2777777777777778
#
# This works for all of the metrics defined in this module. However, because 
# some functions (for example, search functions such as 
# @@hill_climb_stochastic) may take a metric as a parameter, metrics 
# themselves are actually defined as class level variables which are procs.
#
#  > MM.ocm
#  => #<Proc:0x00000100947d38@mm.rb:240 (lambda)>
#
# The "dist_" approach described above is syntactic sugar for the following, 
# more direct call:
#
#  > MM.ocm.call(m,n)
#  => 0.2777777777777778
#
# Though the module provides default settings, each of the metrics is 
# configurable with the use of DistConfig objects. The following example
# turns scaling off (metrics use absolute scaling by default):
#
#  > cfg = MM::DistConfig.new({:scale => :none})
#  => #<MM::DistConfig:0x0000010094b320 @scale=:none, [...] >
#  > MM.dist_ocm(m,n,cfg)
#  => 1.6666666666666667
#
# As mentioned above, a number of the methods in this module are generative. 
# That is, they create morphs which fulfill certain constraints concerning 
# their positions in metric space. For example, to find a point a certain
# distance d away from another point v1 given a metric dist_func:
#
#  > MM.find_point_at_distance({:v1 => m, :d => 0.5, :dist_func => MM.ocm})
#  => NArrayfloat4: [ 1.5, 4.0, 1.0, 3.5 ] 
#  > MM.dist_ocm(m,_)
#  => 0.5
#
# Note that as of this writing, this does not always work. Tweaking calls 
# to ensure their plausibility and double-checking the results are both still
# necessary, as these functions do not return warnings when their results
# are sub-optimal:
#
#  MM.find_point_at_distance({:v1 => m, :d => 0.8, :dist_func => MM.ocm})
#  => NArrayfloat4: [ 4.225, 1.275, 4.275, 1.325 ]
#  > MM.dist_ocm(m,_)
#  => 0.5499999999999999
#
#
#  Some current problems:
#    - Right now MM works with Numeric values only. We want it to work with anything.
#    - Related issue: Contour is currently based on the - operator when it should based 
#      on the <=> operator.
#    - The ordered 2-combinations function only operates on one dimension.
#    - Would like to have a nice function for generating N-dimensional spaces with 
#      several local minima and maxima for testing relative performance of search 
#      algorithms.
#    - Search functions should live in their own separate library.
#
#      
module MM
  include Math
  require 'narray'
  include NMath
  
  # Configuration object for distance functions
  class DistConfig
    attr_accessor :scale, :order, :inter_delta, :intra_delta, :int_func, :ic_calc, :mod
    
    # Creates a new DistConfig object. Takes a configuration hash with the 
    # following keys:
    #
    # [+:scale+] One of +:absolute+, +:relative+, or +:none+. Sets the scaling
    #            method to be applied.
    # [+:order+] The number of times the intra_delta proc is applied to morphs
    #            before their comparison. The default settings (i.e. 
    #            :intra_delta => :abs_diff, :order => 1) correspond to taking
    #            the first discrete derivative of each morph and comparing 
    #            the results.
    # [+:inter_delta+] The delta proc for determining the difference between
    #                  elements of the two morphs. Default: DELTA_FUNCTIONS[:abs_diff]
    # [+:intra_delta+] The delta proc for determining the difference between
    #                  elements of a single morph. Default: DELTA_FUNCTIONS[:abs_diff]
    # [+:int_func+] The interval function proc. See the INTERVAL_FUNCTIONS
    #               constant docs for details. Default: INTERVAL_FUNCTIONS[:plus_one]
    # [+:ic_calc+] The interval class calculation proc. See the IC_FUNCTIONS
    #              constant docs for details. Default: IC_FUNCTIONS[:mod]
    # [+:mod+] The modulus of the numbering system. Used only in interval class
    #          calculations. Default: 12
    #
    def initialize(opts = {})
      @scale       = opts[:scale]       || :absolute
      @order       = opts[:order]       || 1
      @inter_delta = opts[:inter_delta] || MM::DELTA_FUNCTIONS[:abs_diff]
      @intra_delta = opts[:intra_delta] || MM::DELTA_FUNCTIONS[:abs_diff]
      @int_func    = opts[:int_func]    || MM::INTERVAL_FUNCTIONS[:plus_one]
      @ic_calc     = opts[:ic_calc]     || MM::IC_FUNCTIONS[:mod]
      @mod         = opts[:mod]         || 12
    end
  end

  ##
  # :singleton-method: dist_euclidean
  # Simple euclidean distance.
  #
  @@euclidean = ->(m, n, config = nil) do
    # The config argument is only there so the signature is correct
    # for functions that require it, e.g. angle.
    (((m - n) ** 2).sum) ** 0.5
  end
  
  
  #######################
  # Direction Metrics
  #
  # "Direction metrics measure /contour/ differences between morphs. A morph 
  #  and some 'contour preserving' distortion are the same under direction 
  #  metrics. They are listed here from least to most sensitive." MM p. 311
  #

  ## 
  # :singleton-method: dist_uld
  # Unordered Linear Direction, MM p. 312, 316
  # 
  # "[M]easures the differences in average 'up-ness, 'down-ness and equal-
  # ness between two morphs. [...] The ULD is a statistical comparison of
  # linear interval contour, independent of the /corresponding/ respective
  # intra-morphological orders of two morphs. [...] Under the ULD, morphs
  # which 'go up a lot' (linearly) will be closer to other that 'go up a 
  # lot', even if they do not go up in the same places. ULD values range 
  # from [0,1], with a grain of 1 / ((L - 1) * 2)." MM p. 311-12
  #
  # This is the version with scaling for morphs of possibly unequal length
  # given on MM p. 316.
  # :doc:
  @@uld = ->(m, n, config = nil) do
    if config.nil?
      config = self::DistConfig.new
      config.intra_delta = :-.to_proc
    end

    sgn_m = self.sgn(m, config.order, config.intra_delta, config.int_func)
    sgn_n = self.sgn(n, config.order, config.intra_delta, config.int_func)
    sum = 0
    (-1..1).each do |v|
      sum += ((sgn_m.eq(v).sum.to_f / (m.total - 1)) - (sgn_n.eq(v).sum.to_f / (n.total - 1))).abs
    end
    sum.to_f / 2
  end


  ##
  # :singleton-method: dist_old
  # Ordered Linear Direction, MM p. 312
  #
  # "The OLD measures the percentage of different contour values between 
  # corresponding linear intervals. That is, if delta(M)[1] 'goes up' where
  # delta(N)[1] 'goes down' or stays the same, the sum of the direction
  # dissimilarities is incremented. For morphs with the same linear contour,
  # OLD = 0. If two morphs differ in every place (linearly), OLD = 1" 
  # MM p. 313
  #
  @@old = ->(m, n, config = nil) do
    if config.nil?
      config = self::DistConfig.new
      config.intra_delta = :-.to_proc
    end

    sgn_m = self.sgn(m, config.order, config.intra_delta, config.int_func)
    sgn_n = self.sgn(n, config.order, config.intra_delta, config.int_func)
    
    scale_factor = 1
    scale_factor = (m.total - 1) if config.scale == :absolute

    sgn_m.ne(sgn_n).sum.to_f / scale_factor
  end

  ##
  # :singleton-method: dist_ocd
  # Ordered Combinatorial Direction, MM p. 313.
  #
  # "The OCD is the combinatorial version of the OLD. It is the most 
  # discriminating of the direction metrics. The OCD measures the complete,
  # cell-by-cell network of contour similarity between two morphs. The OCD
  # closely reflects melodic perception, tracking the difference between the
  # combinatorial contour of two melodies." MM p. 314
  #
  # The description of diff() (here NArray.ne()) in MM relies on the 
  # generalized interval function. Here instead we are assuming 
  # that we want all of the 2-combinations of items in the vector.
  # We are also assuming we just want the delta between two items, which
  # is necessarily a 1st order calculation.
  #
  # In the event we want something else, we can create that something else
  # and then take the distance from it, rather than providing an interval 
  # function.
  #
  @@ocd = ->(m, n, config = nil) do
    if config.nil?
      config = self::DistConfig.new
      config.intra_delta = :-.to_proc
    end

    # This converts to normal Ruby arrays, and so will be slow.
    m_combo = self.ordered_2_combinations(m.to_a)
    n_combo = self.ordered_2_combinations(n.to_a)

    m_sgn = NArray.to_na(m_combo.map { |a, b| self.sgn_single(config.intra_delta.call(a,b)) })
    n_sgn = NArray.to_na(n_combo.map { |a, b| self.sgn_single(config.intra_delta.call(a,b)) })

    m_sgn.ne(n_sgn).sum.to_f / m_combo.size
  end

  ##
  # :singleton-method: dist_ucd
  # Unordered Combinatorial Direction, MM p. 314, 316.
  #
  # "The UCD compares the statistics of combinatorial 'up/equal/down-ness' of
  # each morph. It does not discern similarities in /corresponding/ 
  # intervals. In general, OCD >= UCD, since the OCD is more sensitive. [...]
  # UCD values range from [0,1] with a grain of 1 / (L[m] * 2)." MM p. 314-16
  #
  # The equation given in MM is misleadingly similar to that given for the
  # ULD metric, but they are not alike; this one is combinatorial.
  #
  # This is the version given for morphs of possibly unequal length, given on
  # MM p. 316
  #
  @@ucd = ->(m, n, config = nil) do
    if config.nil?
      config = self::DistConfig.new
      config.intra_delta = :-.to_proc
    end

    # This converts to normal Ruby arrays, and so will be slow.
    m_combo = self.ordered_2_combinations(m.to_a)
    n_combo = self.ordered_2_combinations(n.to_a)

    m_sgn_map = m_combo.map { |a, b| self.sgn_single(config.intra_delta.call(a,b)) }
    n_sgn_map = n_combo.map { |a, b| self.sgn_single(config.intra_delta.call(a,b)) }

    m_sgn = NArray.to_na(m_sgn_map)
    n_sgn = NArray.to_na(n_sgn_map)

    sum = 0
    (-1..1).each do |v|
      value = ((m_sgn.eq(v).sum.to_f / m_combo.size) - (n_sgn.eq(v).sum.to_f / n_combo.size)).abs
      sum += value
    end
    sum.to_f / 2
  end
  
  
  ################################################################
  # Magnitude Metrics
  #
  
  #
  # A function which returns procs that are themselves magnitude metrics.
  #
  # This function takes a symbol, either +:combinatorial+ or +:linear+, and a
  # proc containing a description of the metric equation. For usage examples,
  # please see the implementation of the four magnitude metrics included here.
  #
  def self.get_mag_metric(style = :combinatorial, post_proc)
    ->(m, n, config = self::DistConfig.new) {
      if style == :combinatorial
        m_combo = self.ordered_2_combinations(m.to_a)
        n_combo = self.ordered_2_combinations(n.to_a)
        m_diff = NArray.to_na(m_combo.map { |a,b| config.intra_delta.call(a,b) })
        n_diff = NArray.to_na(n_combo.map { |a,b| config.intra_delta.call(a,b) })
        #puts "m_combo: #{m_combo.to_a.to_s}"
        #puts "n_combo: #{n_combo.to_a.to_s}"
      elsif style == :linear
        m_combo, n_combo = nil, nil
        m_diff = self.vector_delta(m, config.order, config.intra_delta, config.int_func)
        n_diff = self.vector_delta(n, config.order, config.intra_delta, config.int_func)
      end
      
      #puts "m_diff: #{m_diff.to_a.to_s}"
      #puts "n_diff: #{n_diff.to_a.to_s}"

      scale_factor, inner_scale_m, inner_scale_n = 1, 1, 1

      if config.scale == :absolute
        the_max = [m_diff.max, n_diff.max].max 
        scale_factor = the_max unless the_max == 0
      elsif config.scale == :relative
        inner_scale_m = m_diff.max unless m_diff.max == 0
        inner_scale_n = n_diff.max unless n_diff.max == 0
      elsif config.scale == :maxint_squared
        root_of_squared_differences = ((m_diff - n_diff)**2)**0.5
        scale_factor = root_of_squared_differences.max unless root_of_squared_differences.max == 0
      end

      post_proc.call(config.intra_delta, config.inter_delta, m_diff, n_diff, m_combo, n_combo, 
                     inner_scale_m, inner_scale_n, scale_factor)
    }
  end
  
  ##
  # :singleton-method: dist_olm
  # Ordered Linear Magnitude, generalized. MM p. 305, 318-319
  #
  # "The OLM measures the average difference between corresponding intervals
  # in two morphs. [...] Unlike direction metrics, magnitude metrics are, by
  # definition, unscaled. Since intervals are not bounded, the OLM in its 
  # unscaled form yields indefinitely lage values." MM p. 319
  #
  #--
  #
  # TODO: Clarify this scaling stuff.
  #
  # Absolute scaling: Multiplies the denominator in the averaging operation by 
  #   the max delta (across both vectors). Puts the final average in terms of
  #   % of this max delta.
  # 
  # Relative scaling: Divides each delta vector by its own max. Puts deltas in 
  #   terms of % of the largest delta in each delta vector. I.e. relative 
  #   scaling considers augmented/diminished scalings of vectors to be 
  #   equivalent. E.g. [0,2,6,12] and [0,1,3,6] have distance 0.
  #
  @@olm = self.get_mag_metric(:linear,
    ->(intra_delta, inter_delta, m_diff, n_diff, m_combo, n_combo, 
            inner_scale_m, inner_scale_n, scale_factor) {
      inter_diff = inter_delta.call(m_diff.to_f / inner_scale_m, n_diff.to_f / inner_scale_n)
      #puts "inter_diff: #{inter_diff.to_a.to_s}"
      inter_diff.sum.to_f / (inter_diff.total * scale_factor).to_f
    }
  )

  ##
  # :singleton-method: dist_ulm
  # Unordered Linear Magnitude, generalized. MM p. 304, 320
  #
  # "[T]he ULM measures the /difference of the average/ intervals of two 
  # morphs, whereas the OLM measures the /average difference of corresponding 
  # intervals/ of two morphs. [...] The ULM is not sensitive to intervallic 
  # order, generating a space in which morphs are 'closer' to each other than
  # in OLM-space. In general, OLM >= ULM. Since length(M) need not equal 
  # length(N), the ULM does not require equal length morphs." MM p. 320-21
  #
  #--
  #
  #  TODO: Clarify scaling.
  #
  @@ulm = self.get_mag_metric(:linear,
    ->(intra_delta, inter_delta, m_diff, n_diff, m_combo, n_combo, 
            inner_scale_m, inner_scale_n, scale_factor) {
      inter_delta.call(m_diff.sum.to_f / (m_diff.total * inner_scale_m), 
                       n_diff.sum.to_f / (n_diff.total * inner_scale_n) ).to_f / scale_factor
    }
  )
  
  ##
  # :singleton-method: dist_ocm
  # Ordered Combinatorial Magnitude, MM p. 323
  #
  # "The OCM is the combinatorial version of the OLM. [...] The OCM measures 
  # the average cell-by-cell difference between two absolute magnitude 
  # matrices of equal length morphs." MM p. 324
  #
  #--
  #
  # TODO: This quote probably isn't enough.
  #
  @@ocm = self.get_mag_metric(:combinatorial,
    ->(intra_delta, inter_delta, m_diff, n_diff, m_combo, n_combo, 
            inner_scale_m, inner_scale_n, scale_factor) {
      sum = inter_delta.call(m_diff.to_f / inner_scale_m, 
                             n_diff.to_f / inner_scale_n).sum
      sum.to_f / (m_combo.size * scale_factor)
    }
  )

  ##
  # :singleton-method: dist_ucm
  # Unordered Combinatorial Magnitude, MM p. 325
  #
  # "The UCM is the combinatorial version of the ULM, the difference between
  # average combinatorial intervals. [...] The UCM is a useful statistical 
  # measure. Like the ULM, it does not require that length(m) = length(n)."
  # MM p. 325
  #
  #--
  #
  # TODO: This quote probably isn't enough.
  #
  @@ucm = self.get_mag_metric(:combinatorial,
    ->(intra_delta, inter_delta, m_diff, n_diff, m_combo, n_combo, 
            inner_scale_m, inner_scale_n, scale_factor) {
    #puts "total m: #{(m_diff.to_f / inner_scale_m).sum}"
    #puts "total n: #{(n_diff.to_f / inner_scale_n).sum}"
    #puts "inner scale m: #{inner_scale_m}"
    #puts "inner scale n: #{inner_scale_n}"
    #puts "avgd m: #{(m_diff.to_f / inner_scale_m).sum / (m_combo.size * scale_factor)}"
    #puts "avgd n: #{(n_diff.to_f / inner_scale_n).sum / (n_combo.size * scale_factor)}"
    inter_delta.call((m_diff.to_f / inner_scale_m).sum / (m_combo.size * scale_factor),
                     (n_diff.to_f / inner_scale_n).sum / (n_combo.size * scale_factor)).abs
    }
  )

  
  ##
  # :singleton-method: dist_ic
  # Interval Class Metric, meta-interval form. MM p. 303
  #
  # This differs from the OLM. It calculates interval classes
  # and scales according to the max interval class given a mod.
  # (So for mod 12, the scaling factor is 6.)
  #
  #--
  #
  # TODO: Can this be written in terms of the OLM above?
  #
  @@ic = ->(m, n, config = self::DistConfig.new) do
    scale_factor = 1
    if scale
      scale_factor = mod / 2
    end

    m_diff = self.vector_delta(m, 1, config.intra_delta)
    n_diff = self.vector_delta(n, 1, config.intra_delta)

    inter_diff = config.inter_delta.call(m_diff, n_diff).collect {|x| config.ic_calc.call(x, mod)} 
    inter_diff.sum / (inter_diff.total * scale_factor).to_f
  end


  #
  # Create a multimetric. _MM_ p. 342.
  #
  # Create a metric which combines other metrics. Each item in the +metrics+ 
  # array should be a hash with 3 keys: +:metric+, +:config+, and +:weight+.
  #
  # If weights are omitted, simple averaging is assumed.If configs are 
  # omitted, the default is used. 
  #
  # If the symbol +:none+ is used as a config, the config argument is omitted 
  # in the call to the metric. This is useful if the metric being used is 
  # itself a multimetric.
  #
  def self.get_multimetric(metrics)
    ->(m, n) {  # A config argument would be meaningless.
      top = 0
      bottom = 0
      metrics.each do |metric_hash| #|metric, config = self::DistConfig.new, weight = 1|
        metric = metric_hash[:metric]
        config = metric_hash[:config] || self::DistConfig.new
        weight = metric_hash[:weight] || 1
        
        top += (config == :none) ? (metric.call(m, n) * weight) : (metric.call(m, n, config) * weight)
        bottom += weight
      end
      top.to_f / bottom.to_f
    }
  end
  
  # @@opposition = ->(m, n, v_origin, dist_func = nil, config = self::DistConfig.new)
  # Creating a metric this way would leave us with the wrong type signature,
  # so it couldn't be used in the same way as other metrics. So instead, we:
  
  #
  # Create an opposition metric.
  #
  # Return a lambda which is a metric of the extent to which two points oppose
  # each other around a specified "origin" or third point. This is the scaled 
  # angle relative to two reference points.
  #
  def self.get_opposition_metric(v_origin = nil, dist_func = nil, config = self::DistConfig.new)
    ->(v_ref, v, conf = config) {
      angle = self.angle(v, v_ref, v_origin, dist_func, conf)
      angle / PI
    }
  end
  
  #
  # Make each distance metric lambda readable from the outside, i.e.
  # attr_reader behavior, and also provide sugar for avoiding .call
  # (use dist_#{metric}() instead).
  #
  [:euclidean, :olm, :ulm, :ic, :uld, :old, :ocd, :ucd, :ocm, :ucm, :mag, :ucm_2,
   :hill_climb, :hill_climb_stochastic].each do |sym|
    class_eval(<<-EOS, __FILE__, __LINE__)
      unless defined? @@#{sym}
        @@#{sym} = nil
      end

      def self.#{sym}
        @@#{sym}
      end

      def #{sym}
        @@#{sym}
      end
      
      def self.dist_#{sym}(*args)
        @@#{sym}.call(*args)
      end
    EOS
  end

  #########################################################################
  # Helper functions
  #

  #
  # Magnitude of change between elements in array m. For meta-interval action, 
  # specify your own delta proc. These take the nth order discrete derivative 
  # of the input vector.
  #
  def self.vector_delta(m, order = 1, delta = nil, int_func = nil)
    # Always go from 0 to length - 2, i.e. stop 1 index short of the end.
    # The interval function should seek forward, rather than starting at
    # index 1 and finishing at the final index.

    # So, to get the 1st order discrete derivative, set delta to :-
    # and provide an interval function which generates an array = m[1...m.length]

    if order < 0
      raise "Order must be >= 0"
    elsif order == 0
      return m
    end

    if (delta.class != Symbol && delta.class != Proc && !delta.nil?) || (int_func.class != Symbol && int_func.class != Proc && !int_func.nil?)
      raise "Delta and int_func must be either Symbol or Proc."
    end
    delta = delta.to_proc if delta.class == Symbol

    delta = MM::DELTA_FUNCTIONS[:abs_diff] if delta.nil?
    int_func = MM::INTERVAL_FUNCTIONS[:plus_one] if int_func.nil?

    compare = int_func.call(m)
    # If the compare vector is shorter than m, we assume
    # it is structure such that the missing element at the 
    # end of m is incorporated somehow into the compare vector.
    # E.g. for the default m[i] - m[i+1] approach described in self.
    res = delta.call(m[0...compare.total], compare)
    if order == 1
      return res
    else
      return self.vector_delta(res, order - 1, delta, int_func)
    end
  end

  #
  # Get the interval class of an interval given a modulus.
  #
  # Avoids use of the % operator with negative numbers for compatibility with 
  # both Ruby and NArray definitions.
  #
  def self.interval_class(interval, mod = 12)
    modded = interval.abs % mod
    modded = mod - modded if modded > mod / 2
    modded
  end

  #
  # Number of possible pairwise intervals in a morph. MM p. 307
  # Corresponds to the second-order binomial coefficient of morph length.
  #
  def self.num_pairs(m)
    (m.total**2 - m.total) / 2
  end

  #
  # The contour function, sgn. MM p. 311
  # The sgn function returns:
  #
  #    1 where m[i] > m[j] ("goes down")
  #    0 where m[i] = m[j] ("stays same")
  #   -1 where m[i] < m[j] ("goes up")
  #
  # The formulation in MM is somewhat confusing because it includes
  # the generalized delta symbol, indicating any notion of difference
  # may be plugged in. It seems, though, that "upness" and "downness"
  # necessarily imply a certain notion of difference, i.e. a delta function
  # different from the standard magnitude (abs) function whose use is
  # implicit in most of self.
  #
  # A clue about how this was originally implemented is in the name of
  # the function "sgn" or "sign" and that /negative/ means up and pos
  # means down. This indicates subtraction was used. In the spirit of the
  # generalized delta function in MM, any delta function may be passed in.
  #
  def self.sgn(m, order = 1, intra_delta = :-.to_proc, int_func = nil)
    # We usually want raw difference here, not the magnitude
    deltas = self.vector_delta(m, order, intra_delta, int_func)
    ((deltas > 0).to_f - (deltas < 0)).to_i
  end

  #
  # A version of the sgn function for a single integer as opposed to a whole 
  # vector. I.e., pass in the output from a delta function (e.g. :-) into
  # this, and it will give you results in concordance with the sgn function
  # described above.
  #
  def self.sgn_single(i)
    return -1 if i < 0
    return 1 if i > 0
    0
  end

  #
  # Ordered 2-combinations
  # Provides all of the 2-combinations of an array
  # in an order corresponding to a "flattened" comparison matrix.
  #
  def self.ordered_2_combinations(n)
    combinations = []
    n.each_index do |i|
      ((i + 1)..(n.size - 1)).each do |j|
        combinations << [n[i], n[j]]
      end
    end
    combinations
  end

  #
  # A hash containing delta functions. These functions are for calculating the
  # difference between adjacent values in a morph. Procs contained in this 
  # hash are used in DistConfig configuration objects, and affect the way 
  # metrics are calculated.
  #
  # The included functions are as follows (headings are hash keys):
  #
  # [+:abs_diff+] The absolute value of the difference between elements.
  # [+:raw_diff+] The raw difference between elements.
  # [+:ratio+] The radio of the first element over the second, as a decimal.
  # [+:squared_difference+] The raw difference squared.
  # [+:root_of_squared_difference+] The root of the raw difference squared.
  # [+:huron+] The consonance value of the interval between elements. Uses
  #            David Huron's aggregate dyadic consonance measure.
  #
  # For more information on delta functions, see _MM_ p. 310.
  #
  DELTA_FUNCTIONS = Hash[
    :abs_diff => lambda { |a,b| (a - b).abs },
    :raw_diff => :-.to_proc,
    :ratio => lambda { |a, b| a / b.to_f },
    :squared_difference => lambda { |a, b| (a - b)**2 },
    :root_of_squared_difference => lambda { |a, b| ((a - b)**2)**0.5 },
    :huron => ->(a, b) {
      #              U     m2/M7   M2/m7  m3/M6  M3/m6  P4/P5   A4/d5
      huronTable = [ 0,   -1.428, -0.582, 0.594, 0.386, 1.240, -0.453 ]
      output = (a.to_i - b.to_i).abs
      return huronTable[MM.interval_class(output, 12)] if output.class == Fixnum
      output.collect { |x| huronTable[MM.interval_class(x, 12)] }
    }
  ]
  
  #
  # A hash containing interval functions. Interval functions are procs which
  # take a vector and return a corresponding vector where values with the same
  # index in both vectors are understood as adjacent for the purposes of 
  # ordered metrics. See _MM_ p. 304 for more information.
  #
  # The two functions included here do not come close to exhausting the 
  # possibilities described in _MM_. In particular, at present there is no
  # implemention of the idea of a variable _adjacency_ _interval_.
  #
  # The included functions are as follows (headings are hash keys):
  #
  # [+:plus_one+] The 1st discrete derivative.
  # [+:mean+] Compares each value to the mean of all values.
  #
  INTERVAL_FUNCTIONS = Hash[
    :plus_one => lambda do |m|
      m[1...m.total]          # Default: as in the 1st discrete derivative; use all but the 1st element
    end,
    :mean => lambda do |m|
      NArray.float(m.total).fill!(m.mean)
    end
  ]
  
  #
  # A hash containing interval class functions. Procs contained herein are for 
  # determining the _class_ of a given interval, as opposed to its value
  # in absolute terms. For example, the interval between C and G is always
  # a perfect 5th, regardless of what octave the C and the G are in. The
  # absolute intervals between [60, 67] and [60, 79] are different (7 vs. 19), 
  # but their interval class is the same (5).
  #
  # The included functions are as follows (headings are hash keys):
  #
  # [+:mod+] Uses a modulo system, with a default modulus of 12, corresponding
  #          to the equal-tempered chromatic scale.
  #
  IC_FUNCTIONS = Hash[
    :mod => lambda do |interval, mod = 12|
      self.interval_class(interval, mod)   # Default: Get the mod 12 interval class.
    end
  ]

  #
  # Provides the angle between two vectors with respect to the origin or 
  # another vector (v3).
  #
  def self.angle_euclidean(v1, v2, v3 = nil)
    if !v3.nil?
      v1 -= v3
      v2 -= v3
    end
    Math.acos(dot(v1, v2) / (length(v1) * length(v2)))
  end

  #
  # Get the angle between two vectors given a distance function, e.g. one of 
  # the metrics described above.
  #
  # This only makes sense for certain kinds of metrics. At a minimum, the 
  # metric should satisfy the triangle inequality (i.e. it should be a 
  # metric); beyond that, there are other requirements which I do not yet 
  # entirely understand.
  #
  # Question: Should scaling be turned off here by default?
  # Additional question: What does this even mean?
  #
  def self.angle(v1, v2, v3 = nil, dist_func = nil, config = self::DistConfig.new)
    v3 = NArray.int(v1.total).fill!(0) if v3.nil?
    return 0.0 if (v1 == v3 || v2 == v3)      # not sure if this is strictly kosher
    a = dist_func.call(v1, v3, config)
    b = dist_func.call(v2, v3, config)
    c = dist_func.call(v1, v2, config)
    cos_c = (a**2 + b**2 - c**2).to_f / (2 * a * b)
    Math.acos(cos_c)
  end

  # The Euclidean norm of the vector. Its magnitude.
  #--
  # TODO: This probably needs a less confusing name.
  def self.length(v)
    dot(v,v)**0.5
  end

  # The dot product of two vectors.
  def self.dot(v1, v2)
    (v1 * v2).sum
  end

  # Convert degrees to radians.
  def self.deg2rad(d)
    (d * PI) / 180
  end

  # Convert radians to degrees.
  def self.rad2deg(r)
    (180 * r).to_f / PI
  end

  # Simple linear interpolation in Euclidean space
  def self.interpolate(v1, v2, percent = 0.5)
    ((v2.to_f - v1.to_f) * percent.to_f) + v1.to_f
  end
  
  # Get an array of points along a continuous interpolation between v1 and v2
  # with a given sampling interval, specified as a decimal percentage value.
  # So with interval = 0.1, there'd be 11 samples. (0 and 1 are included.)
  def self.interpolate_path(v1, v2, interval)
    raise "Interval must be > 0 and < 1." if interval < 0.0 || interval > 1.0
    percent = 0.0
    path = []
    while percent < 1.0
      path << self.interpolate(v1, v2, percent)
      percent += interval
    end
    path
  end
  
  # Get an array of points along a continuous interpolation between v1 and v2
  # with a sampling interval determined by the total desired number of points.
  # This method will always include the starting and ending points.
  def self.interpolate_steps(v1, v2, num_steps)
    raise "Number of steps must be > 1" if num_steps <= 0.0
    return [v1, v2] if num_steps == 2
    
    interval = 1.0 / (num_steps - 1)
    current_percent = interval
    path = [v1]
    
    (num_steps - 2).times do |step|
      path << self.interpolate(v1, v2, current_percent)
      current_percent += interval
    end
    
    path << v2
    path
  end
  
  # Upsample a vector, adding members using linear interpolation.
  def self.upsample(v1, new_size)
    return v1 if v1.size == new_size
    raise "Upsample can't downsample" if new_size < v1.size

    samples_to_insert = new_size - v1.size
    possible_insertion_indexes = v1.size - 2
    samples_per_insertion_index = Rational(samples_to_insert, possible_insertion_indexes + 1.0)
    
    count      = Rational(0,1)
    prev_count = Rational(0,1)
    hits       = 0
    new_vector = []
    
    0.upto(possible_insertion_indexes) do |i|
      count += samples_per_insertion_index

      int_boundaries_crossed = count.floor - prev_count.floor
      
      if int_boundaries_crossed >= 1
        hits += int_boundaries_crossed
        # next line leaves off the last step of the interpolation
        # because it will be added during the next loop-through
        new_vector.concat(self.interpolate_steps(v1[i], v1[i + 1], 2 + int_boundaries_crossed)[0..-2])
      else
        new_vector << v1[i]
      end
      prev_count = count
    end
    NArray.to_na(new_vector << v1[-1])    # add the last member (which is not an insertion point)
  end

  #######################################################################
  #
  # Generate vectors based on other vectors
  #

  # Uniformly random point on the euclidean unit n-sphere, scale by radius r
  # Approach from: http://mathworld.wolfram.com/HyperspherePointPicking.html
  def self.sphere_rand_euclidean(n, r = 1.0)
    deviates = NArray.float(n).randomn!
    scale = r / ((deviates**2).sum)**0.5
    scale * deviates
  end
  

  #
  # Find v2 a given distance d from v1 using a given metric and search algorithm.
  # Possibly grossly inefficient. Might hang forever.
  #
  def self.find_point_at_distance(opts)
    v1          = opts[:v1]
    d           = opts[:d]
    dist_func   = opts[:dist_func]
    config      = opts[:config]      || self::DistConfig.new
    search_func = opts[:search_func] || @@hill_climb_stochastic
    search_opts = opts[:search_opts] || {}
    allow_duplicates = opts[:allow_duplicates]
    
    if config == :none
      climb_func = ->(test_point) {
        (dist_func.call(v1, test_point) - d).abs
      }
    else
      climb_func = ->(test_point) {
        (dist_func.call(v1, test_point, config) - d).abs
      }
    end
    search_opts[:climb_func] = climb_func
    search_opts[:start_vector] = v1
    begin
      newpoint = search_func.call(search_opts)
    end while (newpoint.to_a.uniq != newpoint.to_a) && !allow_duplicates
    newpoint
  end
  
  #
  # Find a collection of points a given distance from a vector.
  # 
  # Points on the surface of an n-sphere, but these points will 
  # not necessarily be uniformly distributed.
  #
  def self.set_at_distance(opts)
    set_size     = opts[:set_size]     || 10
    max_failures = opts[:max_failures] || 1000
        
    set = []
    failures = 0
    while set.size < set_size && failures < max_failures
      candidate = self.find_point_at_distance(opts)
      if set.include?(candidate)
        failures += 1
      else
        set << candidate
      end
    end
    set
  end
  
  #
  # Create a path of mutants from one vector to another.
  #
  # This path should compromise between directness in the metric
  # space and directness in Euclidean pitch space. 
  #
  # Results may differ dramatically depending on how the distance 
  # metric is scaled! Adjust euclidean_tightness to taste. High tightness
  # values deemphasize metric space in favor of euclidean space. Typical
  # values would be between 0 and 1.
  #
  # When using low tightness values, it may make sense to cheat at the end.
  #
  # Setting euclidean_tightness to 0 and turning cheating off produces
  # interesting results; perfect in metric space, but wildly divergent in 
  # euclidean space.
  #
  # Some other possibilities: make euclidean_tightness a function or array
  # E.g. it could be tight at the beginning and end but stray in the middle
  #
  def self.metric_path(opts)
    v1     = opts[:v1]
    v2     = opts[:v2]
    metric = opts[:metric]
    config = opts[:config] || MM::DistConfig.new
    steps  = opts[:steps]  || 10
    cheat  = opts[:cheat]  || false
    euclidean_tightness = opts[:euclidean_tightness] || 1.0
    allow_duplicates    = opts[:allow_duplicates]    || true
    search_func = opts[:search_func] || MM.hill_climb_stochastic
    search_opts = opts[:search_opts] || {}
    print_stats = opts[:print_stats] || false
    
    total_distance = (config == :none) ? metric.call(v1, v2) : metric.call(v1, v2, config)
    total_euclidean_distance = MM.euclidean.call(v1,v2)
    inc = total_distance.to_f / steps
    
    puts "total_distance: #{total_distance} inc: #{inc}"
    
    path = [v1]
    steps.times do |step|
      i = step + 1
      climb_func = ->(test_point) {
        dist_v1 = (config == :none) ? metric.call(v1, test_point) : metric.call(v1, test_point, config)
        dist_v2 = (config == :none) ? metric.call(v2, test_point) : metric.call(v2, test_point, config)
        target_v1 = inc * i
        target_v2 = total_distance - (inc * i)
        (target_v1 - dist_v1).abs + (target_v2 - dist_v2).abs + (euclidean_tightness * MM.euclidean.call(test_point, MM.interpolate(v1, v2, i.to_f / steps)))
      }
      #puts "\nNew hill climb, starting from #{path[-1].to_a.to_s}"
      search_opts[:climb_func] = climb_func
      search_opts[:start_vector] = path[-1]

      begin
        newpoint = search_func.call(search_opts)
      end while (newpoint.to_a.uniq != newpoint.to_a) && !allow_duplicates
      path << newpoint

      if print_stats && (vectors_at_each_step == 1)
        puts " - newpoint: #{newpoint.to_a.to_s}"
        dist_v1a = (config == :none) ? metric.call(v1, newpoint) : metric.call(v1, newpoint, config)
        dist_v2a = (config == :none) ? metric.call(v2, newpoint) : metric.call(v2, newpoint, config)
        target_v1a = inc * i
        target_v2a = total_distance - (inc * i)
        off_v1a = (target_v1a - dist_v1a).abs
        off_v2a = (target_v2a - dist_v2a).abs
        off_euc = euclidean_tightness * MM.euclidean.call(newpoint, MM.interpolate(v1, v2, i.to_f / steps))
        puts "target_v1: #{target_v1a} target_v2: #{target_v2a}"
        puts "d_v1: #{metric.call(v1, newpoint, config)} d_v2: #{metric.call(v2, newpoint, config)}"
        puts "off_v1: #{off_v1a} off_v2: #{off_v2a} off_euc: #{off_euc}"
        puts "Eval: #{climb_func.call(newpoint)}"
      end
    end
    path[-1] = v2 if cheat    # Because the hill climb doesn't always make it all the way...
    path
  end
  

  # Point on an n-sphere of radius r
  # and angle less than a from vector v

  # Generate a set of vectors with a given mean & sd

  # Given a mean and a vector count, generate a vector or set of vectors which
  # causes the mean to become a given value
  def self.vector_for_desired_mean(current_mean, count, desired_mean)
    total = current_mean * count
    desired_total = desired_mean * (count + 1)
    desired_total - total
  end
  
  
  
  
  ##################################################
  # Search functions 
  #
  
  #
  # Given an array of coordinates and a proc, evaluate every coordinate 
  # using the proc. Returns a hash using coordinate vectors as keys and
  # return values from the proc as values.
  #
  def self.exhaustive(opts)
    coords = opts[:coords] # ranges per dimension
    func   = opts[:func]   # function to evaluate each coordinate

    results = {}
    
    coords.each do |coord|
      results[coord] = func.call(coord)
    end
    results
  end
  
  def self.generate_coords(opts)
    ranges = opts[:ranges] # ranges per dimension
    incs   = opts[:incs]   # increments per dimension
    
    raise ArgumentError, "opts[:ranges].size must equal opts[:incs].size" if (ranges.size != incs.size) 
    
    possible_value_matrix = []
    
    ranges.each_index {|i|
      range = ranges[i]
      inc   = incs[i]
      
      arr = []
      range.step(inc) {|x| arr << x}
      possible_value_matrix << arr
    }
    
    unfold_pvm(possible_value_matrix)
  end
  
  # Helper function for self.generate_coordinates above.
  # Unfold the possible value matrix.
  # I.e., get all the possible combinations of coordinates without repetition.
  def self.unfold_pvm(pvm, partial_coord = [])
    coords = []
    if pvm.length > 1
      bottom_dim = pvm[0]
      new_pvm = pvm[1..-1]
      bottom_dim.each {|dim_value|
        new_partial_coord = partial_coord.clone
        new_partial_coord << dim_value
        coords.concat(unfold_pvm(new_pvm, new_partial_coord))
      }
    elsif pvm.length == 1
      new_coords = []
      bottom_dim = pvm[0]
      bottom_dim.each {|dim_value|
        new_coord = partial_coord.clone
        new_coord << dim_value
        new_coords << new_coord
      }
      return new_coords
    end
    coords
  end
  
  #
  # Naive hill climbing approach. 
  # Create vectors which minimize a given funcion.
  #
  @@hill_climb = ->(opts) {
    climb_func       = opts[:climb_func]
    start_vector     = opts[:start_vector]
    epsilon          = opts[:epsilon]          || 0.01
    min_step_size    = opts[:min_step_size]    || 0.1
    start_step_size  = opts[:start_step_size]  || 1.0
    max_iterations   = opts[:max_iterations]   || 1000
    return_full_path = opts[:return_full_path] || false                  

    start_vector = start_vector.to_f
    dimensions = start_vector.total
    step_size = NArray.float(dimensions).fill!(start_step_size)
    candidates = [-1, 0, 1]
    current_point = start_vector.dup
    path = [start_vector]
    max_iterations.times do |iteration|
      #puts "- Iteration #{iteration}"
      current_point_cache = current_point.dup
      current_result = climb_func.call(current_point)
      (0...dimensions).each do |i|
        #puts "Testing dimension #{i}"
        best_candidate = nil
        best_score = nil
        candidates.each do |candidate|
          current_point[i] += step_size[i] * candidate
          test_score = climb_func.call(current_point)
          #puts "Testing: #{current_point.to_a.to_s}\t\tScore: #{test_score}"
          current_point[i] -= step_size[i] * candidate
          if best_score.nil? || test_score < best_score
            best_score = test_score
            best_candidate = candidate
          end
        end
        if best_candidate == 0 && step_size[i] > min_step_size
          # Coarse movement on this dimension isn't helping, so
          # lower the size...
          step_size[i] = step_size[i] * 0.5
          step_size[i] = min_step_size if step_size[i] < min_step_size
          #puts "Lowering step size on dimension #{i} to #{step_size[i]}..."
        else
          current_point[i] += step_size[i] * best_candidate
        end
      end
      path << current_point.dup
      test_score = climb_func.call(current_point)
      if test_score < (0 + epsilon) && test_score > (0 - epsilon)
        #puts "Success at iteration #{iteration}"
        break
      end
      if current_point == current_point_cache
        # We can't move...
        #puts "Aborting climb at iteration #{iteration} :("
        break
      end
    end
    return path if return_full_path
    current_point
  }
  
  ##
  # :singleton-method: hill_climb_stochastic
  #
  # A stochastic hill climbing algorithm. Probably a bad one.
  #
  # Takes a hash with the following keys:
  #
  # [+:climb_func+] A proc. This is the function the algorithm will attempt 
  #                 to minimize.
  # [+:start_vector+] The starting point of the search.
  # [+:epsilon+] Anything point below this distance will be considered a 
  #              successful search result. Default: 0.01
  # [+:min_step_size+] The minimum step size the search algorithm will use. Default: 0.1
  # [+:start_step_size+] The starting step size. Default: 1.0
  # [+:max_iterations+] The number of steps to take before giving up.
  # [+:return_full_path+] If true, this will return every step in the search.
  # [+:step_size_subtract+] If set to a value, every time the algorithm needs
  #                         to decrease the step size, it will decrement it
  #                         by this value. Otherwise it will divide it by 2.0.
  #
  @@hill_climb_stochastic = ->(opts) {
    climb_func         = opts[:climb_func]
    start_vector       = opts[:start_vector]
    epsilon            = opts[:epsilon]            || 0.01
    min_step_size      = opts[:min_step_size]      || 0.1
    start_step_size    = opts[:start_step_size]    || 1.0
    max_iterations     = opts[:max_iterations]     || 1000
    return_full_path   = opts[:return_full_path]   || false
    step_size_subtract = opts[:step_size_subtract]

    start_vector = start_vector.to_f
    dimensions = start_vector.total
    step_size = start_step_size
    candidates = [-1, 0, 1]
    current_point = start_vector.dup
    path = [start_vector]
    max_iterations.times do |iteration|
      #puts "- Iteration #{iteration}"
      current_point_cache = current_point.dup
      current_result = climb_func.call(current_point)
      # Generate a collection of test points with scores: [point, score]
      test_points = []
      (dimensions * 5).times do |j| # num test points = num dimensions * 3
        # pick a deviation for each dimension
        deviation = NArray.float(dimensions)
        dimensions.times { |d| deviation[d] = candidates.sample * step_size }
        new_point = current_point + deviation
        test_points << [new_point, climb_func.call(new_point)]
      end
      
      test_points.sort_by! { |x| x[1] }
      winner = test_points[0]
      if winner[1] < current_result
        current_point = winner[0] 
        path << current_point.dup
      end
      
      test_score = climb_func.call(current_point)
      if test_score < (0 + epsilon) && test_score > (0 - epsilon)
        #puts "Success at iteration #{iteration}"
        break
      end
      if current_point == current_point_cache && step_size > min_step_size
        # We didn't get any good results, so lower the step size
        if step_size_subtract
          step_size -= step_size_subtract 
        else
          step_size = step_size * 0.5
        end
        step_size = min_step_size if step_size < min_step_size
        #puts "Lower step size to #{step_size}"
      elsif current_point == current_point_cache && step_size <= min_step_size
        # We didn't get any good results, and we can't lower the step size
        #puts "Aborting climb at iteration #{iteration} :("
        break
      end
    end
    return path if return_full_path
    current_point
  }
  
  
end
