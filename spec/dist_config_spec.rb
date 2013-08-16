require_relative '../mm.rb'

require 'minitest/autorun'

describe MM::DistConfig do
  describe "with no options passed" do
    let(:dist_config) { MM::DistConfig.new }

    it "should have default scaling" do
      dist_config.scale.must_equal :absolute
    end
    
    it "should have default order" do
      dist_config.order.must_equal 1
    end
    
    it "should have default mod" do
      dist_config.mod.must_equal 12
    end
    
    # TODO: This tests no edge cases
    it "should have default delta functions" do
      dist_config.inter_delta.call(5, 2).must_equal 3
      dist_config.inter_delta.call(2, 5).must_equal 3
      
      dist_config.intra_delta.call(5, 2).must_equal 3
      dist_config.intra_delta.call(2, 5).must_equal 3
    end
    
    it "should have default interval function" do
      dist_config.int_func.call(NArray[2,3,4,5]).must_equal NArray[3,4,5]
    end
    
    it "should have default interval class calculation function" do
      dist_config.ic_calc.call(6).must_equal 6
      dist_config.ic_calc.call(18).must_equal 6
    end
  end
end
