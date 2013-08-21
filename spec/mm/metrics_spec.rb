require 'spec_helper'
require 'minitest/autorun'
require 'mm'

describe MM do
  let(:m) {NArray[1,5,12,2,9,6]} # first point
  let(:n) {NArray[7,6,4,9,8,1]} # second point
  let(:o) {NArray[5,3,6,1,4,3]} # third point for testing triangle inequality
  let(:dist_config) {nil}
  
  describe "direction metrics" do
    describe "OLD Metric" do
      let(:subject) {MM.old}
      let(:expected) {0.8}
      it_behaves_like "a metric"
    end
    
    describe "OCD Metric" do
      let(:subject) {MM.ocd}
      let(:expected) {0.666}
      it_behaves_like "a metric"
    end
    
    describe "ULD Metric" do
      let(:subject) {MM.uld}
      let(:expected) {0.4}
      it_behaves_like "a metric"
    end
    
    describe "UCD Metric" do
      let(:subject) {MM.ucd}
      let(:expected) {0.2666}
      it_behaves_like "a metric"
    end
  end
  
  describe "magnitude metrics" do
    describe "OLM Metric" do
      let(:subject) {MM.olm}
      describe "with scale: none" do
        let(:expected) {4.6}
        let(:dist_config) {MM::DistConfig.new :scale => :none}
        it_behaves_like "a metric"
      end
      
      describe "with scale: absolute" do
        let(:expected) {0.46}
        let(:dist_config) {MM::DistConfig.new :scale => :absolute}
        it_behaves_like "a metric"
      end
      
      describe "with scale: relative" do
        let(:expected) {0.44}
        let(:dist_config) {MM::DistConfig.new :scale => :relative}
        it_behaves_like "a metric"
      end
    end
    
    describe "OCM Metric" do
      let(:subject) {MM.ocm}
      
      describe "with scale: none" do
        let(:expected) {3.6}
        let(:dist_config) {MM::DistConfig.new :scale => :none}
        
        it_behaves_like "a metric"
      end
      
      describe "with scale: absolute" do
        let(:expected) {0.327}
        let(:dist_config) {MM::DistConfig.new :scale => :absolute}
        
        it_behaves_like "a metric"
      end
      
      describe "with scale: relative" do
        let(:expected) {0.366}
        let(:dist_config) {MM::DistConfig.new :scale => :relative}
        
        it_behaves_like "a metric"
      end
    end
    
    describe "ULM Metric" do
      let(:subject) {MM.ulm}
      
      describe "with scale: none" do
        let(:expected) {3.0}
        let(:dist_config) {MM::DistConfig.new :scale => :none}
        
        it_behaves_like "a metric"
      end
      
      describe "with scale: absolute" do
        let(:expected) {0.3}
        let(:dist_config) {MM::DistConfig.new :scale => :absolute}
        
        it_behaves_like "a metric"
      end
      
      describe "with scale: relative" do
        let(:expected) {0.1628}
        let(:dist_config) {MM::DistConfig.new :scale => :relative}
        
        it_behaves_like "a metric"
      end
    end
    
    describe "UCM Metric" do
      let(:subject) {MM.ucm}
      
      describe "with scale: none" do
        let(:expected) {1.6}
        let(:dist_config) {MM::DistConfig.new :scale => :none}
        
        it_behaves_like "a metric"
      end
      
      describe "with scale: absolute" do
        let(:expected) {0.14545}
        let(:dist_config) {MM::DistConfig.new :scale => :absolute}
        
        it_behaves_like "a metric"
      end
      
      describe "with scale: relative" do
        let(:expected) {0.025}
        let(:dist_config) {MM::DistConfig.new :scale => :relative}
        
        it_behaves_like "a metric"
      end
    end
  end
  
  describe "euclidean metric" do
    let(:subject) {MM.euclidean}
    let(:expected) {13.27}
    
    it_behaves_like "a metric"
  end
  
  describe "ic metric" do
    let(:subject) {MM.ic}
    # I'm not sure this is the correct result
    let(:expected) {0.7666}
    
    it_behaves_like "a metric"
  end
end