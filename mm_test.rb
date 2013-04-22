require './mm.rb'
require 'test/unit'

#
# These tests are for the metrics themselves...
#
class MMTest < Test::Unit::TestCase
  
  def test_olm
    # p. 331-334
    m = NArray[1,5,12,2,9,6]
    n = NArray[7,6,4,9,8,1]
    assert_in_delta(4.6,  MM.olm.call(m,n, MM::DistConfig.new(:scale => :none))    , 0.01)
    assert_in_delta(0.46, MM.olm.call(m,n, MM::DistConfig.new(:scale => :absolute)), 0.01)
    assert_in_delta(0.44, MM.olm.call(m,n, MM::DistConfig.new(:scale => :relative)), 0.01)
  end
  
  def test_old
    # p. 331-334
    m = NArray[1,5,12,2,9,6]
    n = NArray[7,6,4,9,8,1]
     # MM says this value should be 0.2. But I think it should actually be 0.8...
    assert_in_delta(0.8, MM.old.call(m,n), 0.01)
    o = NArray[5,3,6,1,4]
    p = NArray[3,6,1,4,2]
    assert_in_delta(1.0, MM.old.call(o,p), 0.01)
  end
  
  def test_ocm
    # p. 323-324
    o = NArray[1,6,2,5,11]
    p = NArray[3,15,13,2,9]
    assert_in_delta(5.2, MM.ocm.call(o,p, MM::DistConfig.new(:scale => :none)), 0.01)
    
    # p. 331-334
    m = NArray[1,5,12,2,9,6]
    n = NArray[7,6,4,9,8,1]
    assert_in_delta(3.6,   MM.ocm.call(m,n, MM::DistConfig.new(:scale => :none))    , 0.01)
    assert_in_delta(0.327, MM.ocm.call(m,n, MM::DistConfig.new(:scale => :absolute)), 0.01)
    assert_in_delta(0.366, MM.ocm.call(m,n, MM::DistConfig.new(:scale => :relative)), 0.01)  
  end
  
  def test_ocd
    # p. 331-334    
    m = NArray[1,5,12,2,9,6]
    n = NArray[7,6,4,9,8,1]
    assert_in_delta(0.666, MM.ocd.call(m,n), 0.001)
  end
  
  def test_ulm
    # Old comments suggest there is something funny about these, from p. 330
    m = NArray[30,35,45,43]
    n = NArray[1,25,3,9]
    o = NArray[30,31,28,33]
    #assert_in_delta(4.67, MM.ulm(m,n), 0.01)
    #assert_in_delta(1.56, MM.ulm(m,o), 0.01)
    #assert_in_delta(6.23, MM.ulm(n,o), 0.01)

    # p. 331-334    
    m = NArray[1,5,12,2,9,6]
    n = NArray[7,6,4,9,8,1]
    assert_in_delta(3.0,    MM.ulm.call(m,n, MM::DistConfig.new(:scale => :none))    , 0.01)  
    assert_in_delta(0.3,    MM.ulm.call(m,n, MM::DistConfig.new(:scale => :absolute)), 0.01)  
    assert_in_delta(0.1628, MM.ulm.call(m,n, MM::DistConfig.new(:scale => :relative)), 0.0001)
  end
  
  def test_ucm
    # p. 331-334  
    m = NArray[1,5,12,2,9,6]
    n = NArray[7,6,4,9,8,1]
    assert_in_delta(1.6,     MM.ucm.call(m,n, MM::DistConfig.new(:scale => :none))    , 0.01)
    assert_in_delta(0.14545, MM.ucm.call(m,n, MM::DistConfig.new(:scale => :absolute)), 0.00001)
    assert_in_delta(0.025,   MM.ucm.call(m,n, MM::DistConfig.new(:scale => :relative)), 0.001)
  end
  
  def test_uld
    # p. 331-334  
    m = NArray[1,5,12,2,9,6]
    n = NArray[7,6,4,9,8,1]
    assert_in_delta(0.4, MM.uld.call(m,n), 0.01)
    
    o = NArray[5,3,6,1,4]
    p = NArray[3,6,1,4,2]
    assert_in_delta(0.0, MM.uld.call(o,p), 0.01)
  end
  
  def test_ucd
    # p. 331-334  
    m = NArray[1,5,12,2,9,6]
    n = NArray[7,6,4,9,8,1]
    assert_in_delta(0.2666, MM.ucd.call(m,n), 0.0001)    
  end
  
  ###################
  # Nesting Metrics 
  #
  # Ensure that metrics can be used as delta functions for other metrics
  # TODO: Write tests to every combination of metric and intra/inter delta
  # 
  
  # MM.dist_ocm
  # - using the UCD as intra-delta
  # - using abs_diff as inter_delta
  def test_ucd_as_intra_delta
    # Drawn from MM, p. 315
    q = NArray[5,3,7,6]
    r = NArray[2,1,2,1]
    s = NArray[8,3,5,4]
    # Setting up two different orderings of the intervals
    m = NArray[q, r, s]
    n = NArray[r, s, q]
    # Setting up the OCD intra_delta with hard-coded scaling
    ocd_proc = ->(a, b, config = MM::DistConfig.new(:scale => :relative)) do
      MM.dist_ocd(a, b, config)
    end
    # Expected results:
    # Using the OCD as intra_delta
    # m_combos = [0.5, 0.33, 0.33]
    # n_combos = [0.33, 0.5, 0.33]
    # 
    # (scaled over 0.5)
    # OCM diffs = [0, 0.3333, 0.3333]
    assert_in_delta(0.222, MM.dist_ocm(m, n, MM::DistConfig.new(:intra_delta => ocd_proc, :scale => :none)), 0.001)
    assert_in_delta(0.666, MM.dist_ocm(m, n, MM::DistConfig.new(:intra_delta => ocd_proc, :scale => :absolute)), 0.001)
  end
  
  def test_exhaustive
    opts_gen = {
      :ranges => [(0..2), (0..2), (0..2)],
      :incs   => [1, 1, 1]
    }
    
    coords = MM.generate_coords(opts_gen)       
           
    opts_search = {
      :coords => coords, 
      :func   => ->(coord) {
        coord == [1, 1, 1]
      }
    }
            
    res = {[0, 0, 0]=>false, [0, 0, 1]=>false, [0, 0, 2]=>false, [0, 1, 0]=>false, [0, 1, 1]=>false, [0, 1, 2]=>false, [0, 2, 0]=>false, [0, 2, 1]=>false, [0, 2, 2]=>false, [1, 0, 0]=>false, [1, 0, 1]=>false, [1, 0, 2]=>false, [1, 1, 0]=>false, [1, 1, 1]=>true, [1, 1, 2]=>false, [1, 2, 0]=>false, [1, 2, 1]=>false, [1, 2, 2]=>false, [2, 0, 0]=>false, [2, 0, 1]=>false, [2, 0, 2]=>false, [2, 1, 0]=>false, [2, 1, 1]=>false, [2, 1, 2]=>false, [2, 2, 0]=>false, [2, 2, 1]=>false, [2, 2, 2]=>false}
    
    assert_equal(res, MM.exhaustive(opts_search))
  end
  
end

#
# These tests are for helper functions, etc...
#
class MMHelperTest < Test::Unit::TestCase
  def test_unfold_pvm
    pvm0 = [[0, 1, 2], [0, 1, 2], [0, 1, 2]]
    
    res0 = [[0, 0, 0], [0, 0, 1], [0, 0, 2],
            [0, 1, 0], [0, 1, 1], [0, 1, 2],
            [0, 2, 0], [0, 2, 1], [0, 2, 2],
                                           
            [1, 0, 0], [1, 0, 1], [1, 0, 2],
            [1, 1, 0], [1, 1, 1], [1, 1, 2],
            [1, 2, 0], [1, 2, 1], [1, 2, 2],
                                           
            [2, 0, 0], [2, 0, 1], [2, 0, 2],
            [2, 1, 0], [2, 1, 1], [2, 1, 2],
            [2, 2, 0], [2, 2, 1], [2, 2, 2]]
           
    pvm1 = [[0,1]]
    res1 = [[0],[1]]
    
    pvm2 = [[0,1],[0,1]]
    res2 = [[0,0], [0,1], [1,0], [1,1]]
           
    assert_equal(res0, MM.unfold_pvm(pvm0))
    assert_equal(res1, MM.unfold_pvm(pvm1))
    assert_equal(res2, MM.unfold_pvm(pvm2))
  end
end
