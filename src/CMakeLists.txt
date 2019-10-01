FILE(GLOB holstein_src "*.cc" "*.cpp" "external/*.cc")

list(REMOVE_ITEM holstein_src "._*.cc")  # I'm not convinced that this works...
list(REMOVE_ITEM holstein_src "._*.cpp")  # I'm not convinced that this works...

ADD_LIBRARY(holstein_lib ${holstein_src})
TARGET_LINK_LIBRARIES(holstein_lib ${ROOT_LIBRARIES})