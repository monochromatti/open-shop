import json,pathlib,hashlib,sys,collections
stages={
 'screen':('57cf943a6f1ee0f6f251166faede10968208ebeb','37729035609',7,1),
 'combinations':('3eb507390e6f14408715995d072679d222c13a0a','37730112283',5,1),
 'refinement':('c4397c6909516c9ebd26f368323296adcf357f4d','37732208118',5,1),
 'confirmation':('c4397c6909516c9ebd26f368323296adcf357f4d','37733769105',3,2),
}
measurements=[];checks=[]
for stage,folder in [item.split('=',1) for item in sys.argv[1:]]:
 revision,workflow,profiles,repeats=stages[stage]
 paths=list(pathlib.Path(folder).rglob('summary.json'))
 for path in sorted(paths):
  rows=json.loads(path.read_text())
  assert len(rows)==profiles*repeats,(stage,path,len(rows))
  for key in ('case_sha256','seed_controls_sha256','manifest_sha256','source_sha256','experiment_sha256'):
   assert len({r[key] for r in rows})==1,(stage,path,key)
  lower=max(r['feasible_lower_bound'] for r in rows)
  for row in rows:
   assert not row.get('experiment_error'),(stage,path,row.get('experiment_error'))
   assert row['accepted'] and row['start_audit']['valid']
   assert row['global_bound']>=lower-1e-6
   assert all(p['upper']>=lower-1e-6 for p in row['global_bound_trajectory'])
   assert not row['bound_event_errors']
   sep=row.get('targeted_supports')
   assert sep is None or (not sep['errors'] and not sep['infeasible_flags'])
   native=row.pop('scip_statistics')
   compact={k:native.get(k) for k in ('status','timing','origprob','presolvedprob','tree','root','solution','lp')}
   compact['propagators']={k:v for k,v in native.get('propagator',{}).get('plugins',{}).items() if k in ('obbt','genvbounds','nonlinear')}
   compact['nonlinear_handlers']={k:v for k,v in native.get('nlhdlr',{}).get('plugins',{}).items() if k in ('bilinear','default','perspective')}
   row['native_statistics']=compact
   stats=path.parent/(row['profile']+'-'+row['commitment']+'-'+str(row['repeat'])+'.log.statistics.json')
  row['native_statistics_sha256']=hashlib.sha256(stats.read_bytes()).hexdigest()
   row.update(stage=stage,revision=revision,workflow='https://github.com/monochromatti/open-shop/actions/runs/'+workflow,artifact=path.parent.name)
   measurements.append(row)
  checks.append(dict(stage=stage,artifact=path.parent.name,records=len(rows),best_audited_objective=lower))
assert measurements, 'no measurements found'
assert {r['manifest_sha256'] for r in measurements}=={'d0079266b836f68bbf9d5fc159d4f571d3bf3b6ef104dc2c8c3448bacfbbfa97'}
assert all(r['bound_event_monotonic'] and r['threads']==1 and r['blas_threads']==1 for r in measurements)
print(json.dumps(dict(schema_version=1,baseline_revision='7743702366ad96ce2ecfa24881c66a44a887c9a5',validation=checks,measurements=measurements),indent=2))
