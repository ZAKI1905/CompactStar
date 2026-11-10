"""Generate an isolated, instrumented Stepper for read-only arithmetic evidence.

Only class naming and passive Apply begin/end callbacks differ from production.
This class is never used by a production trajectory or linked into its library.
"""
def generate(root,build):
 trace=(root/'tests/diagnostics/checkpoint_replay.cpp').read_text()
 trace=trace[trace.index('struct TraceContext'):trace.index('std::string Bytes')]
 trace=trace.replace('ActualCheckpointContext context;','std::size_t rows=0;ActualCheckpointContext context;')
 trace=trace.replace('int cell){rhs<<','int cell){++rows;rhs<<')
 assert '++rows' in trace
 trace+=r'''
static TraceContext* audit_context=nullptr;
static std::ostream* audit_stream=nullptr;
static double audit_start_time=0;
static State audit_start_state{};
static std::size_t audit_start_row=0,audit_apply_id=0;
void BindAudit(TraceContext& context,std::ostream& stream){
 audit_context=&context;audit_stream=&stream;audit_apply_id=0;
 stream<<std::setprecision(17)<<"context\tapply\tleft\tright\trow_begin\trow_end\trejected\tx0\te0\tmu0\tx1\te1\tmu1\n";
}
void AuditApplyBegin(double t,const State& y){
 require(audit_context&&audit_stream,"unbound arithmetic audit");
 audit_start_time=t;audit_start_state=y;audit_start_row=audit_context->rows;
}
void AuditApplyEnd(double t,const State& y,std::size_t rejected){
 *audit_stream<<audit_context->id<<'\t'<<++audit_apply_id<<'\t'<<audit_start_time<<'\t'<<t<<'\t'<<audit_start_row<<'\t'<<audit_context->rows<<'\t'<<rejected<<'\t'<<values(audit_start_state)<<'\t'<<values(y)<<'\n';
 audit_stream->flush();require(bool(*audit_stream),"arithmetic boundary write failed");
}
'''
 (build/'trace_context.inc').write_text(trace)
 src=(root/'CompactStar/Physics/BNV/BreakpointStepper.hpp').read_text()
 src=src.replace('BreakpointStepper','AuditStepper')
 old='const int rc=gsl_odeiv2_evolve_apply('
 assert src.count(old)==1
 src=src.replace(old,'::AuditApplyBegin(t,y);\n        '+old)
 old='stats_.rejected+=evolve_->failed_steps-rejected;'
 assert src.count(old)==1
 src=src.replace(old,old+'\n        ::AuditApplyEnd(t,y,evolve_->failed_steps-rejected);')
 (build/'AuditStepper.hpp').write_text(src)
