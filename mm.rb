module MM
  include Math
  require 'narray'
  include NMath
  
  #
  # Configuration objects for distance functions
  #
  class DistConfig
    attr_accessor :scale, :order, :inter_delta, :intra_delta, :int_func, :ic_calc, :mod

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

  #
  # These distance functions expect vector input (I'm using NArray)
  # i.e. a Morph with a single parameter
  #

  #
  # Euclidean distance
  #
  @@euclidean = ->(m, n, config = nil) do
    # The config argument is only there so the signature is correct
    # for functions that require it, e.g. angle.
    (((m - n) ** 2).sum) ** 0.5
  end

  #
  # Absolute magnitude metric, MM p. 299
  #
  @@amm = ->(m, n, config = self::DistConfig.new) do
    result = (m - n).abs.sum 
    result = result / m.total.to_f if scale == :absolute
    result
  end
  
  ################################################################
  # Magnitude Metrics
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
  
  #
  # Ordered Linear Magnitude, generalized. MM p. 305, 318-319
  #
  # The average difference between the deltas.
  #
  # Absolute scaling: Multiplies the denominator in the averaging operation by 
  #   the max delta (across both vectors). Puts the final average in terms of
  #   % of this max delta.
  # 
  # Relative scaling: Divides each delta vector by its own max. Puts deltas in 
  #   terms of % of the largest delta in each delta vector. I.e. relative scaling
  #   considers augmented/diminished scalings of vectors to be equivalent. E.g. 
  #   [0,2,6,12] and [0,1,3,6] have distance 0.
  #
  @@olm = self.get_mag_metric(:linear,
    ->(intra_delta, inter_delta, m_diff, n_diff, m_combo, n_combo, 
            inner_scale_m, inner_scale_n, scale_factor) {
      inter_diff = inter_delta.call(m_diff.to_f / inner_scale_m, n_diff.to_f / inner_scale_n)
      #puts "inter_diff: #{inter_diff.to_a.to_s}"
      inter_diff.sum.to_f / (inter_diff.total * scale_factor).to_f
    }
  )

  #
  # Unordered Linear Magnitude, generalized. MM p. 304, 320
  #
  # The difference between average deltas for each vector.
  # 
  @@ulm = self.get_mag_metric(:linear,
    ->(intra_delta, inter_delta, m_diff, n_diff, m_combo, n_combo, 
            inner_scale_m, inner_scale_n, scale_factor) {
      inter_delta.call(m_diff.sum.to_f / (m_diff.total * inner_scale_m), 
                       n_diff.sum.to_f / (n_diff.total * inner_scale_n) ).to_f / scale_factor
    }
  )
  
  #
  # Ordered Combinatorial Magnitude, MM p. 323
  #
  @@ocm = self.get_mag_metric(:combinatorial,
    ->(intra_delta, inter_delta, m_diff, n_diff, m_combo, n_combo, 
            inner_scale_m, inner_scale_n, scale_factor) {
      sum = inter_delta.call(m_diff.to_f / inner_scale_m, 
                             n_diff.to_f / inner_scale_n).sum
      sum.to_f / (m_combo.size * scale_factor)
    }
  )

  #
  # Unordered Combinatorial Magnitude, MM p. 325
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


  #######################
  # Direction Metrics
  #

  #
  # Unordered Linear Direction, MM p. 312, 316
  # 
  # This is the version with scaling for morphs of possibly unequal length
  # given on MM p. 316.
  #
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


  #
  # Ordered Linear Direction, MM p. 312
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

  #
  # Ordered Combinatorial Direction, MM p. 313.
  #
  #   5 9 3 2       2 5 6 6
  # 5   - + +     2   - - -
  # 9     + +     5     - -
  # 3       +     6       0
  # 2             6
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

  #
  # Unordered Combinatorial Direction, MM p. 314, 316.
  #
  # The equation given in MM is misleadingly similar to that given for the
  # ULD metric, but they are not alike; this one is combinatorial.
  #
  # The situation with passing in a delta function is the same as for OCD
  # above, for the same reasons.
  #
  # This is the version given for morphs of possibly unequal length, given on MM p. 316
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
  
  #
  # Interval Class Metric, meta-interval form. MM p. 303
  #
  # This differs from the OLM. It calculates interval classes
  # and scales according to the max interval class given a mod.
  # (So for mod 12, the scaling factor is 6.)
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
  # Create a multimetric. MM p. 342.
  #
  # Create a metric which combines other metrics.
  # Each item in the metrics array should be a hash
  # with 3 keys: :metric, :config, :weight
  #
  # If weights are omitted, simple averaging is assumed. 
  # If configs are omitted, the default is used. 
  #
  # If the symbol :none is used as a config, the config argument is omitted in the 
  # call to the metric. (E.g. if the metric being used is itself a multimetric.)
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
  
  # Return a lambda which is a metric of the extent to which two points oppose
  # each other around a specified "origin" or third point.
  # This is the scaled angle relative to two reference points...  
  def self.get_opposition_metric(v_origin = nil, dist_func = nil, config = self::DistConfig.new)
    ->(v_ref, v, conf = config) {
      angle = self.angle(v, v_ref, v_origin, dist_func, conf)
      angle / PI
    }
  end
  
  ###
  #
  # Make each distance metric lambda readable from the outside, i.e.
  # attr_reader behavior, and also provide sugar for avoiding .call
  # (use dist_#{metric}() instead).
  #
  [:euclidean, :amm, :olm, :ulm, :ic, :uld, :old, :ocd, :ucd, :ocm, :ucm, :mag, :ucm_2,
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
  # Magnitude of change between elements in array m
  # For meta-interval action, specify your own delta lambda
  # These take the nth order discrete derivative of the input vector
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
  # Get the interval class of an interval given a modulus
  # It's all gotta be Integers.
  #
  # Avoids use of the % operator with negative numbers for
  # compatibility with both Ruby and NArray definitions.
  #
  # For mod = 12, 0-6 should return 0-6 and...
  # 7 -> 5, 8 -> 4, 9 -> 3, 10 -> 2, 11 -> 1, 12 -> 0
  # 13 -> 1, 14 -> 2, ... 19 ->5, 20-> 4
  def self.interval_class(interval, mod = 12)
    modded = interval.abs % mod
    modded = mod - modded if modded > mod / 2
    modded
  end

  #
  # Number of possible pairwise intervals in a morph, MM p. 307
  # Second-order binomial coefficient of morph length.
  #
  def self.num_pairs(m)
    (m.total**2 - m.total) / 2
  end

  #
  # The contour function, sgn. MM p. 311
  # The sgn function returns:
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
  # A version of the sgn function for a single integer
  # as opposed to a whole vector.
  # I.e., pass in the output from a delta function (e.g. :-) into
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
  # Alternative delta functions. MM p. 310.
  #
  DELTA_FUNCTIONS = Hash[
    :abs_diff => lambda { |a,b| (a - b).abs },         # Default: abs of difference
    :raw_diff => :-.to_proc,
    :ratio => lambda { |a, b| a / b.to_f },
    :squared_difference => lambda { |a, b| (a - b)**2 },
    :root_of_squared_difference => lambda { |a, b| ((a - b)**2)**0.5 },
    :huron => ->(a, b) {
      #              U       m2/M7   M2/m7  m3/M6  M3/m6  P4/P5   A4/d5
      huronTable = [ 0,   -1.428, -0.582, 0.594, 0.386, 1.240, -0.453 ]
      output = (a.to_i - b.to_i).abs
      return huronTable[MM.interval_class(output, 12)] if output.class == Fixnum
      output.collect { |x| huronTable[MM.interval_class(x, 12)] }
    }
  ]

  INTERVAL_FUNCTIONS = Hash[
    :plus_one => lambda do |m|
      m[1...m.total]          # Default: as in the 1st discrete derivative; use all but the 1st element
    end,
    :mean => lambda do |m|
      NArray.float(m.total).fill!(m.mean)
    end
  ]
  
  IC_FUNCTIONS = Hash[
    :mod => lambda do |interval, mod = 12|
      self.interval_class(interval, mod)   # Default: Get the mod 12 interval class.
    end
  ]

  # Provides the angle between two vectors with respect to the origin
  # or another vector (v3)
  def self.angle_euclidean(v1, v2, v3 = nil)
    if !v3.nil?
      v1 -= v3
      v2 -= v3
    end
    Math.acos(dot(v1, v2) / (length(v1) * length(v2)))
  end

  # Get the angle between two vectors given a distance function
  # e.g. one of the metrics described above.
  # (This is why the triangle equality is important.)
  #
  # Question: Should scaling be turned off here by default?
  # Additional question: What does this even mean?
  def self.angle(v1, v2, v3 = nil, dist_func = nil, config = self::DistConfig.new)
    v3 = NArray.int(v1.total).fill!(0) if v3.nil?
    return 0.0 if (v1 == v3 || v2 == v3)      # not sure if this is strictly kosher
    a = dist_func.call(v1, v3, config)
    b = dist_func.call(v2, v3, config)
    c = dist_func.call(v1, v2, config)
    cos_c = (a**2 + b**2 - c**2).to_f / (2 * a * b)
    Math.acos(cos_c)
  end

  def self.length(v)
    dot(v,v)**0.5
  end

  def self.dot(v1, v2)
    (v1 * v2).sum
  end

  def self.deg2rad(d)
    (d * PI) / 180
  end

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

  # Find v2 a given distance d from v1 using a given metric and search algorithm.
  def self.find_point_at_distance(v1, d, dist_func, config = self::DistConfig.new, search_func = @@hill_climb_stochastic)
    if config.nil?
      climb_func = ->(test_point) {
        (dist_func.call(v1, test_point) - d).abs
      }
    else
      climb_func = ->(test_point) {
        (dist_func.call(v1, test_point, config) - d).abs
      }
    end
    search_func.call(climb_func, v1)
  end
  
  #
  # Find a collection of points a given distance from a vector.
  # 
  # Points on the surface of an n-sphere, but these points will 
  # not necessarily be uniformly distributed. Also grossly inefficient.
  #
  def self.set_at_distance(opts = {})
    v1           = opts[:v1]
    d            = opts[:d]
    dist_func    = opts[:dist_func]
    config       = opts[:config]       || self::DistConfig.new
    
    search_func  = opts[:search_func]  || @@hill_climb_stochastic
    search_epsilon    = opts[:search_epsilon]    || 0.01
    search_min_step   = opts[:search_min_step]   || 1.0
    search_start_step = opts[:search_start_step] || 4.0
    set_size     = opts[:set_size]     || 10
    max_failures = opts[:max_failures] || 1000
    
    if config == :none
      climb_func = ->(test_point) {
        (dist_func.call(v1, test_point) - d).abs
      }
    else
      climb_func = ->(test_point) {
        (dist_func.call(v1, test_point, config) - d).abs
      }
    end
    
    set = []
    failures = 0
    while set.size < set_size && failures < max_failures
      candidate = search_func.call(climb_func, v1, search_epsilon, search_min_step, search_start_step)
      if !set.include?(candidate)
        set << candidate
      else
        failures += 1
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
  def self.metric_path(v1, v2, metric, config = nil, steps = 10, euclidean_tightness = 1.0, cheat = false, allow_duplicates = true)
    total_distance = config.nil? ? metric.call(v1, v2) : metric.call(v1, v2, config)
    total_euclidean_distance = MM.euclidean.call(v1,v2)
    inc = total_distance.to_f / steps
    
    puts "total_distance: #{total_distance} inc: #{inc}"
    
    path = [v1]
    steps.times do |step|
      i = step + 1
      climb_func = ->(test_point) {
        dist_v1 = metric.call(v1, test_point, config)
        dist_v2 = metric.call(v2, test_point, config)
        target_v1 = inc * i
        target_v2 = total_distance - (inc * i)
        (target_v1 - dist_v1).abs + (target_v2 - dist_v2).abs + (euclidean_tightness * MM.euclidean.call(test_point, MM.interpolate(v1, v2, i.to_f / steps)))
      }
      #puts "\nNew hill climb, starting from #{path[-1].to_a.to_s}"
      
      begin
        newpoint = MM.hill_climb_stochastic.call(climb_func, path[-1], 0.001, 1.0, 4.0)
      end while (newpoint.to_a.uniq != newpoint.to_a) && !allow_duplicates
      puts " - newpoint: #{newpoint.to_a.to_s}"
      dist_v1a = metric.call(v1, newpoint, config)
      dist_v2a = metric.call(v2, newpoint, config)
      target_v1a = inc * i
      target_v2a = total_distance - (inc * i)
      off_v1a = (target_v1a - dist_v1a).abs
      off_v2a = (target_v2a - dist_v2a).abs
      off_euc = euclidean_tightness * MM.euclidean.call(newpoint, MM.interpolate(v1, v2, i.to_f / steps))
      puts "target_v1: #{target_v1a} target_v2: #{target_v2a}"
      puts "d_v1: #{metric.call(v1, newpoint, config)} d_v2: #{metric.call(v2, newpoint, config)}"
      puts "off_v1: #{off_v1a} off_v2: #{off_v2a} off_euc: #{off_euc}"
      puts "Eval: #{climb_func.call(newpoint)}"
      path << newpoint
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
  # Naive hill climbing approach. 
  # Create vectors which minimize a given funcion.
  #
  @@hill_climb = ->(climb_func, start_vector, epsilon = 0.01, min_step_size = 0.1, 
                    start_step_size = 1.0, max_iterations = 1000, return_full_path = false) {
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
  
  #
  # Stochastic hill climbing
  #
  @@hill_climb_stochastic = ->(climb_func, start_vector, epsilon = 0.01, min_step_size = 0.1, 
                      start_step_size = 1.0, max_iterations = 1000, return_full_path = false) {
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
        step_size = step_size * 0.5
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
