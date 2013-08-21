module MM
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
end