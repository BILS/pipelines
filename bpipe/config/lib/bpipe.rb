# A module to manage Bpipe pipelines for the BILS pipeline code base
module Bpipe
  
  # == DESCRIPTION
  # A Pipeline object holds information on all installed Bpipe pipelines.
  class Pipeline
  
    attr_accessor :module_folder , :pipeline_folder, :pipelines, :modules, :db, :root_folder
  
    # Creates a new pipeline object. Requires BPIPE_LIB to be set 
    def initialize    
      abort "BPIPE_LIB not set, aborting!" unless ENV['BPIPE_LIB']
      
      @root_folder = ENV['BPIPE_LIB'].gsub(/\/modules/, '')
      @module_folder = "#{@root_folder}/modules"
      @pipeline_folder = "#{@root_folder}/pipelines"

      # This order is not random!
      
      @db = []
      _parse_db_file

      @modules = []
      _parse_pipeline_modules     
      
      @pipelines = []
      _get_pipeline_files
    
    end
  
    # Lists all available pipelines
    def listall
      puts "BILS pipelines:"

      self.pipelines.sort_by{|p| p.name }.each do |pipeline|
        puts "----------------"
        puts "Name:\t\t#{pipeline.name}"
        puts "Title:\t\t#{pipeline.title}"
	puts "Input(s):\t#{pipeline.inputs.join('\n')}"
      end    
      exit
    end
	
    # Lists all requirements from all modules
    def list_requirements     
      puts @modules.collect{|m| m.requirements}.flatten.uniq.join("\n")
      exit
    end
    
    # Copies the bpipe config template to the current working directory
    def bpipe_config
      system("cp #{root_folder}/templates/bpipe.config.template .")
      exit
    end

    # Checks the internal database of variables and descriptions to make sure
    # that all variables used by any of the stages/modules is present. Reports
    # any variables that are not yet included. 
    def db_verify
      issues = []
      requirements = self.modules.collect{|m| m.requirements}.flatten.uniq
      requirements.each do |req|
        # If the required stage is not known to the module db, make a noise
	unless self.db.find{|db| db.variable == req }
          warn "The following variable is not defined in config/config_db: #{req}"
          issues << req
        end
      end      
      warn "A total of #{issues.length} problems were found, consider updating the database!" unless issues.empty?
      exit
    end

    private 
	
    # Internal method that generates Pipeline::PipelineFile objects
    # These correspond to workflow logics (in pipelines/)
    def _get_pipeline_files
      Dir.entries(self.pipeline_folder).select{ |file| file.include?(".bpipe") }.each do |file|
        this_pipeline = Bpipe::Pipeline::PipelineFile.new(self.pipeline_folder + "/" + file)
	# The is probably a neater way of doing this, but...:
	this_pipeline.add_modules(@modules)
        this_pipeline.requirements.each do |r|
          db_entry = self.db.find{|entry| entry.variable == r}
          this_pipeline.db << db_entry unless db_entry.nil?
        end
        @pipelines << this_pipeline
      end
    end
  
    # Internal method to compile all modules that are available in this
    # code base.
    def _parse_pipeline_modules
      Dir.entries(self.module_folder).select{|file| file.include?(".groovy") }.each do |file|
        mods = Bpipe::Pipeline::ModuleParser.new(self.module_folder + "/" + file).parse
        mods.each do |mod|
          @modules << mod
        end
      end
      @modules = self.modules.sort_by{|m| m.name }
    end

    # Method to parse a database of pipeline variables and descriptions (config/config_db)
    def _parse_db_file
      
      IO.readlines(@root_folder + "/config/config_db").each do |line|
        elements = line.strip.split("\t")
	if elements.length == 3 
          @db << Bpipe::Pipeline::PipelineVariable.new(elements.shift,elements.shift,elements.shift)
	elsif elements.length == 2
          @db << Bpipe::Pipeline::PipelineVariable.new(elements.shift,elements.shift)
	else
          abort "Entry in config_db has unrecognized structure:\n#{line}"
	end
      end      

    end

    # == DESCRIPTION
    # The ModuleParser class can read module files (modules/*.groovy) and
    # extract bpipe modules from it. Reports Pipeline::PipelineModule objects.
    class ModuleParser < Bpipe::Pipeline
	
      attr_accessor :path, :parsed_modules

      def initialize(file)
        @path = file
        @parsed_modules = []
        _parse_module_file(file)
      end

      def parse
        return @parsed_modules
      end

      private 
	   
      def _parse_module_file(file)
        lines = IO.readlines(file)
        stream = []
        lines.each do |line|
          if line.include?("= {") or line.include?("= segment {")# start of a module
  		    stream = [line.strip]
		  elsif line.match(/^\}$/) # end of a module
  		    stream << line.strip
                    this_module = Pipeline::PipelineModule.new(stream)
                    this_module.filename = file # We store the name of the file holiding this module
                    self.parsed_modules << this_module
  	      else
  		    stream << line.strip
  	      end
        end
      end

    end

    # == DESCRIPTION
    # The PipelineModule class parses individual Bpipe modules and
    # stores information on name and defined requirements for each 
    # module. 
    # Expects the module as an array of lines, which are normally
    # passed from ModuleParser::_parse_module_file
    class PipelineModule < Bpipe::Pipeline
      
      attr_accessor :name , :requirements, :path, :description, :default_path, :filename

      def initialize(lines)
        @name = nil
        @requirements = []
        @description = ""
        @default_path = ""
        @exec = nil
        @filename = nil
        _parse_string(lines)
      end
       
      private
	  
      # Parses the module and stores identified values in designated variables. 
      def _parse_string(lines)
        self.requirements = lines.select{|line| line.include?("requires") and line.include?(":") }.collect{|line| line.split(":")[0].gsub(/^.*requires/, '').strip }
        self.name = lines[0].split("=")[0].strip
      end

    end

    # == DESCRIPTION
    # The PipelineFile class is a parser class that extracts information
    # on individual stages and their requirements from a bpipe pipeline file
    # (pipelines/*.bpipe)
    class PipelineFile < Bpipe::Pipeline

      attr_accessor :title, :name, :commands, :modules, :requirements, :db, :inputs

      def initialize(file)
        @title = nil
        @name = file.split("/")[-1].gsub(/\.bpipe/, '')
        @commands = []
        @db = []
	@inputs = []
        @modules = []
        _parse_pipeline_file(file)
                		
      end

      # Adds modules to this pipeline (from an array of modules, vetted against stages)
      def add_modules(modules)
        self.commands.each do |command|
          this_module = modules.find{|m| m.name == command }
          abort "Pipeline #{self.name} contains a stage name without a matching module: #{command}" unless this_module
          self.modules << this_module unless self.modules.include?(this_module)
        end
      end

       # Adds an entry to the database
      def add_db_entry(entry)
	self.db << entry unless self.db.include?(entry)
      end

      # Provides information about a pipeline
      def describe
        puts "Pipeline information for #{self.name}"
        puts "Description: #{self.title}"
        puts
        puts "Modules used:"
        self.modules.each do |mod|
          puts "\t#{mod.name} (in: #{mod.filename})"
        end
        puts
        puts "Variables required:"
        self.requirements.each do |entry|
          puts "\t#{entry}"
        end
        puts ""
      end
       
      # Writes a config files for this pipeline,
      # written to pipeline.config.template
      # If the DB holds information on a variable,
      # it will be added as comment.
      def write_config

        f = File.new("pipeline.config.template","w+")
        f.puts "// Config file for #{self.name}"
        self.requirements.each do |r|
          db_entry = self.db.find{ |entry| entry.variable == r }
          f.puts
          if db_entry
            f.puts "//#{db_entry.description}"
            f.puts "#{db_entry.variable}=\"#{db_entry.path}\""
          else
            f.puts "//No description available yet"
            f.puts "#{r}=\"\""
          end
        end
        f.close

      end
	  
      def requirements
        return @modules.collect{|m| m.requirements}.flatten.uniq
      end
 
      private
	  
      # Parses a pipeline file and returns information on the stages being
      # used as well as the title
      def _parse_pipeline_file(file)

        lines = IO.readlines(file)

	# Parse the title of the pipeline
        self.title = lines.find{|line| line.include?("title:") }.slice(/\".*\"/)

	# Parse the expected inputs to this pipeline
	self.inputs = lines.select{|l| l.strip.match(/^inputs.*/)}.collect{|l| l.gsub(/\"/, '').split(" ")[1..-1].join(" ") }

	# Clean some of the obious junk
	run_string = lines.collect{|l| l.strip}.join(" ").slice(/\{.*\}/).gsub(/[\[\]\{\}]/, '').gsub(/,/, ' ')

	# Clean decorational elements and extract the functions
	run_string = run_string.gsub(/[\+\*]/, '').gsub(/\"(.*?)\"/, '').gsub(/\.using\((.*?)\)/, '')

	# Store the functions used by this pipeline
	if run_string
        	run_string.split(" ").collect{|e| self.commands << e }
	end
		
      end
    
      # Writes a config file which holds all identified
      # requirements based on the used modules. 
  
    end

    # == DESCRIPTION
    # A class to hold information on variables used
    # across different modules
    class PipelineVariable < Bpipe::Pipeline

      attr_accessor :variable, :description, :path

      def initialize(variable,description,path="")
        @variable = variable
        @description = description
        @path = path
      end

    end
  
  end

end
