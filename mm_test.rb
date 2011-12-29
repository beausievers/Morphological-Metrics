require './mm.rb'
require 'test/unit'

class MMTest < Test::Unit::TestCase
  def test_ocd
    m = NArray[5,9,3,2]
    n = NArray[2,5,6,6]
    o = NArray[5,3,6,1,4]
    p = NArray[3,6,1,4,2]
  end
end

#
# Tests on p. 330:
#
# TODO: What's up with these...?
#

m = NArray[30,35,45,43]
n = NArray[1,25,3,9]
o = NArray[30,31,28,33]

MM.ulm(m,n) # 4.67
MM.ulm(m,o) # 1.56
MM.ulm(n,o) # 6.23


# Tests on p. 331 - 334

m = NArray[1,5,12,2,9,6]
n = NArray[7,6,4,9,8,1]

MM.olm.call(m,n, MM::DistConfig.new(false))      # 4.6
MM.olm.call(m,n, MM::DistConfig.new("absolute")) # 0.46
MM.olm.call(m,n, MM::DistConfig.new("relative")) # 0.44

MM.ocm.call(m,n, MM::DistConfig.new(false))      # 3.6
MM.ocm.call(m,n, MM::DistConfig.new("absolute")) # 0.327
MM.ocm.call(m,n, MM::DistConfig.new("relative")) # 0.366
o = NArray[1,6,2,5,11]
p = NArray[3,15,13,2,9]
MM.ocm.call(o,p, MM::DistConfig.new(false))      # 5.2    from p. 323-324

MM.ulm.call(m,n, MM::DistConfig.new(false))      # 3.0
MM.ulm.call(m,n, MM::DistConfig.new("absolute")) # 0.3
MM.ulm.call(m,n, MM::DistConfig.new("relative")) # 0.1628

MM.ucm.call(m,n, MM::DistConfig.new(false))      # 1.6
MM.ucm.call(m,n, MM::DistConfig.new("absolute")) # 0.14545
MM.ucm.call(m,n, MM::DistConfig.new("relative")) # 0.025

MM.uld.call(m,n) # 0.4

MM.old.call(m,n) # 0.2  or so says MM. But I think this value should actually be 0.8...
o = NArray[5,3,6,1,4]
p = NArray[3,6,1,4,2]
MM.old.call(o,p) # 1.0
MM.uld.call(o,p) # 0.0

MM.ucd.call(m,n) # 0.2666
MM.ocd.call(m,n) # 0.6666