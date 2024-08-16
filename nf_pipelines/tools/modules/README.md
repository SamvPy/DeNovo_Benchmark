# De novo engine modules

Each file should be a stand-alone process which takes as input:
- an MGF channel 
- a configuration channel (optional)
- a serializer file

Inside the process, the tool is run and an output is generated:
- result channel with naming convention <filename>.<tool>.<any-extension>
- a log channel if required
- the serializer file

IMPORTANT! The file should be named as <filename>.<tool>.<any-extension>

This is important for compatibility with the next processes in the pipeline as the filename will be used,
to group files together.

Recommended:
The storeDir directive is recommended, so that the pipeline can be run partially, without requiring the long runtime
of a redundant computation, while providing a directory to backtrace the original outputs in case something goes wrong during parsing.

If this behaviour is not required, publishDir can be used. The naming of the file is recommended as <filename_mgf>.<any-extension>
