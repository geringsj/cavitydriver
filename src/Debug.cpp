/**
 * I hate myself for doing this, but we need to do this if we want to set
 * the output of the debug and log defines from the parseArguments function
 * in IO. 
 * TODO: Find a nicer way to do this!!!
 */
bool m_debug;
bool m_log_err;
bool m_log_warn;
bool m_log_info;