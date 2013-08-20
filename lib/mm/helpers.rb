module MM
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
end