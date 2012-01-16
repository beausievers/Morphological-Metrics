require './mm.rb'
require 'test/unit'

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
end

