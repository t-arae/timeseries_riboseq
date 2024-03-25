module Fs

import Base.Filesystem; const fs = Filesystem

dir_create(path) = fs.mkpath(path)
dir_ls(path = fs.pwd(); full::Bool = false) = fs.readdir(path, join = full)

path_ext(path) = fs.splitext(path)[2]
# fs_dir_ls() |> x -> fs_path_ext.(x)
path_ext_rm(path) = fs.splitext(path)[1]
# fs_dir_ls() |> x -> fs_path_ext_rm.(x)
path_ext_set(path, ext) = path * ext
path_ext_rm_set(path, ext) = path_ext_rm(path) * ext
# fs_path_ext_rm_set("hoge.sam", ".sort.bam")

path_file(path) = fs.basename(path)
# fs_dir_ls() |> x -> fs_path_file.(x)
path_dir(path) = fs.dirname(path)
# fs_dir_ls(full = true) |> x -> fs_path_dir.(x)

path_filter(path, regex) = occursin(regex, path) ? path : nothing

end # module