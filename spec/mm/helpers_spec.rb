require 'spec_helper'
require 'minitest/autorun'
require 'mm'

describe MM do
  let(:m) { NArray[2,5,6,1] }
  
  describe "::vector_delta" do
    describe "with default arguments" do
      it "returns absolute difference" do
        MM.vector_delta(m).must_equal NArray[3,1,5]
      end
      
      it "returns same data type" do
        MM.vector_delta(m.to_f).must_equal NArray[3.0,1.0,5.0]
      end
    end
    
    describe "with higher order" do
      it "returns higher order" do
        MM.vector_delta(m, 2).must_equal NArray[2,4]
        MM.vector_delta(m, 3).must_equal NArray[2]
      end
      
      it "raises error if order is too high" do
        ->{MM.vector_delta(m, 4)}.must_raise IndexError
      end
    end
    
    describe "with alternative delta" do
      let(:delta) { MM::DELTA_FUNCTIONS[:raw_diff] }
      
      it "uses the delta passed" do
        MM.vector_delta(m, 1, delta).must_equal NArray[-3, -1, 5]
      end
    end
  end
  
  describe "::interval_class" do
    it "finds mod-12 ics" do
      MM.interval_class(18).must_equal 6
      MM.interval_class(-3).must_equal 3
      MM.interval_class(-12).must_equal 0
      MM.interval_class(9).must_equal 3
    end
    
    it "finds alternative mod ics" do
      MM.interval_class(18, 17).must_equal 1
      MM.interval_class(6, 3).must_equal 0
      MM.interval_class(8, 3).must_equal 1
    end
  end
  
  describe "::num_pairs" do
    it "finds number of possible pairs in a morph" do
      MM.num_pairs(m).must_equal 6
    end
  end
  
  describe "::sgn" do
    describe "with default arguments" do
      let(:expected) {NArray[-1, -1, 1]}
      
      it "returns contour of a morph" do
        MM.sgn(m).must_equal expected
      end
      
      it "returns int no matter the data type" do
        MM.sgn(m.to_f).must_equal expected
      end
      
      it "includes 0 as a possibility" do
        MM.sgn(NArray[-1, -1, 1]).must_equal NArray[0, -1]
      end
    end
    
    describe "with higher order" do
      # My expectation was that the higher-order sgn function would take the 
      # contour of the contour. It should be noted exactly what's going on, so
      # that there's no confusion.
      it "returns higher order" do
        MM.sgn(m, 2).must_equal NArray[-1, -1]
      end
    end
  end
  
  describe "::sgn_single" do
    describe "when given a single number" do
      it "returns the sign of the number" do
        MM.sgn_single(-100).must_equal -1
        MM.sgn_single(3.5).must_equal 1
        MM.sgn_single(0).must_equal 0
      end
    end
  end
  
  describe "::ordered_2_combinations" do
    describe "when given a one-dimensional vector" do
      let(:vector) { [1,2,3] }
      let(:expected) { [[1,2],[1,3],[2,3]] }
      it "returns all combinations" do
        MM.ordered_2_combinations(vector).must_equal expected
      end
    end
    
    describe "when given a higher-dimensional vector" do
      let(:vector) { [[2,1],[3,2],[4,3]] }
      let(:expected) { [[[2,1], [3,2]], [[2,1], [4,3]], [[3,2], [4,3]]] }
      it "returns combinations of the outermost element" do
        MM.ordered_2_combinations(vector).must_equal expected
      end
    end
    
    describe "when given non-numeric values" do
      let(:vector) { ["hi there", "how are you?", "fine, thanks"] }
      let(:expected) { [["hi there", "how are you?"], ["hi there", "fine, thanks"], ["how are you?", "fine, thanks"]]}
      it "returns combinations of the elements" do
        MM.ordered_2_combinations(vector).must_equal expected
      end
    end
  end
  
  describe "DELTA_FUNCTIONS" do
    describe ":abs_diff" do
      it "returns the absolute difference" do
        MM::DELTA_FUNCTIONS[:abs_diff].call(-5, 3).must_equal 8
        MM::DELTA_FUNCTIONS[:abs_diff].call(3, -5).must_equal 8
      end
    end
    
    describe ":raw_diff" do
      it "returns the raw difference" do
        MM::DELTA_FUNCTIONS[:raw_diff].call(-5, 3).must_equal -8
        MM::DELTA_FUNCTIONS[:raw_diff].call(3, -5).must_equal 8
      end
    end
    
    describe ":ratio" do
      it "returns the float ratio" do
        MM::DELTA_FUNCTIONS[:ratio].call(-5, 3).must_be_close_to -1.6667, 0.001
        MM::DELTA_FUNCTIONS[:ratio].call(3, -5).must_be_close_to -0.6, 0.001
      end
    end
    
    describe ":squared_difference" do
      it "returns the squared difference" do
        MM::DELTA_FUNCTIONS[:squared_difference].call(-5, 3).must_equal 64
        MM::DELTA_FUNCTIONS[:squared_difference].call(3, -5).must_equal 64
      end
    end
    
    describe ":root_of_squared_differences" do
      it "returns the root of squared difference" do
        MM::DELTA_FUNCTIONS[:root_of_squared_difference].call(-5, 3).must_equal 8
        MM::DELTA_FUNCTIONS[:root_of_squared_difference].call(3, -5).must_equal 8
      end
    end
    
    describe ":huron" do
      it "returns the huron table value of the interval" do
        MM::DELTA_FUNCTIONS[:huron].call(-5, 3).must_equal 0.386
        MM::DELTA_FUNCTIONS[:huron].call(3, -5).must_equal 0.386
      end
    end
  end
  
  describe "INTERVAL_FUNCTIONS" do
    describe ":plus_one" do
      let(:vector) { NArray[1, 2, 3] }
      let(:expected) { NArray[2, 3] }
      let(:subject) { MM::INTERVAL_FUNCTIONS[:plus_one] }
      
      it "returns all but the first element" do
        subject.call(vector).must_equal expected
      end
    end
    
    describe ":mean" do
      let(:vector) { NArray[1, 2, 3] }
      let(:expected) { NArray[2, 2, 2] }
      let(:subject) { MM::INTERVAL_FUNCTIONS[:mean] }
      
      it "returns an array of only the mean values" do
        subject.call(vector).must_equal expected
      end
    end
  end
  
  describe "IC_FUNCTIONS" do
    describe ":mod" do
      it "returns the modulo interval class" do
        MM::IC_FUNCTIONS[:mod].call(7).must_equal 5
        MM::IC_FUNCTIONS[:mod].call(-3).must_equal 3
        MM::IC_FUNCTIONS[:mod].call(19).must_equal 5
      end
      
      it "accepts an alternate modulus" do
        MM::IC_FUNCTIONS[:mod].call(19, 15).must_equal 4
      end
    end
  end
end
