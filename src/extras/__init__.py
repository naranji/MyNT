# append parent folder to package path,
# such that custom modules located in parent folder can be found
# NOTE: not really nice, should maybe changed into global package
__path__.append(__path__[0].rpartition('/')[0])