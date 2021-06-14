export exp_begin, exp_cancelled, exp_comment,
	exp_end, exp_filename, exp_prefix,
	exp_reset_timer, exp_save, exp_savefig

include("Dbg.jl")

using Dates, Base.Filesystem, CPUTime, InfoZIP, JLD, Logging, Printf, Plots, Images
using Manopt

function exp_begin(;aux_folders = "", prefix::String = "", postfix::String = "")

global results_folder
global results_prefix
global exp_cancelled_v
global io

exp_cancelled_v = false

dt = now()
runid = Dates.format(dt, "yyyy-mm-dd-HH-MM-SS")
results_folder = "results/" * prefix * runid * postfix

mkpath(results_folder)
results_prefix = results_folder * "/"

io = open(results_prefix * "log.txt", "w+")
logger = SimpleLogger(io)
global_logger(logger)

files = String[]
if aux_folders != ""
    if aux_folders isa String
        aux_folders = [aux_folders]
	elseif aux_folders isa Array{String,1}
    	for i in 1:length(aux_folders)
			files = vcat(files, [aux_folders[i] * "/*.jl"])
    	end
	end
end
files = vcat(files, ["*.jl"])
rp = replace(results_prefix, "/" => "\\") # otherwise error for some reason
#create_zip(rp * "src.zip", files)

exp_reset_timer()

end

function exp_cancelled()
  global exp_cancelled_v
  return exp_cancelled_v
end

function exp_comment(onwhat::String = "")
	global io
	global results_prefix

	if onwhat == ""
		prefix = ""
	else
		prefix = onwhat * ": "
	end

	print("Comment: ")
	msg = readline()
	@printf io "%s\n" "Comment: " * msg

	fid  = open(results_prefix * "_comment.txt", "w+")
	write(fid, prefix * msg * "\n")
	close(fid)

end

function exp_end()
	close(io)
end

function exp_filename(filename::String)
	global results_prefix
  	result = results_prefix * filename
end

function exp_io()
	global io
	return io
end

# function exp_inspect()
# end

function exp_prefix()
 	global results_prefix
  	result = results_prefix
end

function exp_reset_timer()
 	global debug_timer_start
  	debug_timer_start = CPUtic()
end

# function exp_get_time()
# 	global debug_timer_start
# 	time = CPUtime_us()
# 	return (time-debug_timer_start)/1e6
# end

function exp_save(filename::String, varargin...)

	global results_prefix
 	fn = results_prefix * filename * ".jld"

	if !isempty(varargin)
		s = Dict()
		i = 1
		while i + 1 <= length(varargin)
			s[varargin[i]] = varargin[i+1]
			i = i + 2
		end
		save(fn, s)
	else
		# TODO this still does nothing. Q: do we need it?
		# save(fn)
	end
	return fn
end

function exp_savefig(filename::String, save_eps::Bool = false)
	global results_prefix

	savefig(results_prefix * filename * ".png")

	if save_eps
		savefig(results_prefix * filename * ".eps")
	end
end

function exp_savefig(plot_ref::Plots.Plot, filename::String, save_eps::Bool = false)
	global results_prefix

	savefig(plot_ref, results_prefix * filename * ".png")

	if save_eps
		savefig(plot_ref, results_prefix * filename * ".eps")
	end
end

function exp_saveim(img::Array{K,2}, filename::String) where K
	global results_prefix

	save(results_prefix * filename * ".png",img)
end

function exp_savetable(table::String, filename::String) where K
	global results_prefix

	io = open(results_prefix * filename * ".txt","w")
	write(io,table)
	close(io)
end
