module MM
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
end