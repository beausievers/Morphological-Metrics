$: << File.expand_path(File.dirname(__FILE__) + "/../lib")

class Minitest::SharedExamples < Module
  include Minitest::Spec::DSL
  
  def self.find(name)
    @shared_examples ||= {}
    @shared_examples[name]
  end
  
  def self.register(name, &block)
    @shared_examples ||= {}
    @shared_examples[name] = Minitest::SharedExamples.new(&block)
  end
end

class Minitest::Spec
  def self.it_behaves_like(name)
    include(Minitest::SharedExamples.find(name))
  end
  
  def self.shared_examples_for(name, &block)
    ::Minitest::SharedExamples.register(name, &block)
  end
  
  # Implementing shared expectations for all metrics
  shared_examples_for "a metric" do
    
    it "is transitive" do
      subject.call(m, n, dist_config).must_equal subject.call(n, m, dist_config)
    end
    
    it "is the shortest distance between two points" do
      o = NArray[5,3,6,1,4,3]
      subject.call(m, n, dist_config).must_be :<=, (subject.call(m, o, dist_config) + subject.call(o, n, dist_config))
    end
    
    it "returns the correct distance" do
      subject.call(m, n, dist_config).must_be_close_to expected, 0.01
    end
  end
end