import matlab.engine
eng = matlab.engine.start_matlab()
eng.buck_boost_meas(nargout=0) # Expects a file named myMatlabFile.m in the same directory
#eng.quit()


