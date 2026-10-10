from pathlib import Path
import re,json,os,subprocess
base=Path('research_audit/weighted_formal').resolve(); src=Path('formalize_11_12/math/lean/OAI/NumberTheory/DirichletL')
files=['Detector/FinalAssemblyData','Moments/NaturalFixedRaySourceFixedIdeal','Moments/NaturalFixedRaySourceFreeExceptional','Moments/DetectorEnergyInitialState','Moments/DetectorPlainSlotProfile','Moments/DetectorPlainCapacity','Moments/DetectorPlainExceptional','Moments/DetectorPlainMarkedState','Moments/DetectorPlainMomentParameters','ParametersHighDataFine','Detector/DetectorMomentPrimeGeometry','Dictionary/InverseMarkedReferencePrimeGates','Dictionary/InverseMarkedReferenceDeletedEnergy','Dictionary/InverseMarkedReferenceFiberEnergy','Detector/DetectorInverseMarkedField','Detector/DetectorInverseFields','Detector/DetectorPlainUnmarkedField','Detector/DetectorPlainUnmarkedRestrictedField','Detector/DetectorPlainMarkedFineField','Detector/FinalAssemblyMomentInput','Detector/FinalAssemblyCertifiedBands']
contents={}
for f in files:
 s=(src/(f+'.lean')).read_text()
 if f.endswith('FinalAssemblyData'):s=s[:s.index('\ndef SourceData.probe')]+ '\nend SevenEighths.ProbeFinalAssembly\nend\nend OAI\n'
 if f=='ParametersHighDataFine':s=s[:s.index('theorem exists_high_data_fine')]+s[s.index('theorem fine_slot_widths'):]
 if f.endswith('FinalAssemblyMomentInput'):s=s[:s.index('\nlemma polynomial_height_absorption')]+ '\nend SevenEighths.ProbeFinalAssembly\nend\nend OAI\n'
 if f.endswith('FinalAssemblyCertifiedBands'):
  s=s[:s.index('\ndef DetectorCertifiedBands')]+s[s.index('\ntheorem plain_inputs'):s.index('\ntheorem fine_moments_of_certified')]+ '\nend SevenEighths.ProbeFinalAssemblyCertifiedBands\nend\nend OAI\n'
  s=s.replace('detectorMesh (Δ t:ℝ)','detectorMesh (t:ℝ)').replace('{Δ t:ℝ}(hΔ:0<Δ)','{t:ℝ}').replace('detectorMesh Δ','detectorMesh')
 contents[f]=s
names=set()
for s in contents.values():
 for n in re.findall(r'^(?:def|lemma|theorem|structure)\s+([\w.]+)',s,re.M):
  if '.' not in n and n!='SourceData':names.add(n)
replacements={n:'weighted_'+n for n in names}; replacements['SourceData']='WeightedSourceData'
pat=re.compile(r'(?<![\w])('+ '|'.join(sorted(map(re.escape,replacements),key=len,reverse=True))+r')(?![\w])')
record='''import OAI.NumberTheory.DirichletL.ParametersHighData
namespace WeightedQRH
structure MomentData where
  N : ℕ
  ell : Fin N → ℝ
  t : ℝ
  ε : ℝ
  e : ℝ
  t_pos : 0 < t
  t_small : t ≤ 1/100000000
  slots_injective : Function.Injective ell
  slots_bounds : ∀j, 0 < ell j ∧ ell j ≤ 1
  epsilon_pos : 0 < ε
  epsilon_small : ε ≤ 1/1000
  e_pos : 0 < e
end WeightedQRH
'''
(base/'WeightedQRH/MomentData.lean').write_text(record)
prev='MomentData'
for f,s in contents.items():
 target='Moment'+f.split('/')[-1]
 s=pat.sub(lambda m:replacements[m[0]],s)
 s=re.sub(r'(?:Parameters\.)?HighData Δ','WeightedQRH.MomentData',s)
 s=s.replace('3/4+2*Δ','3/4').replace('(3/4 + 2 * Δ)','(3/4)')
 s=s.replace('{Δ : ℝ} ','').replace('{Δ:ℝ}','').replace('variable \n','')
 if f.endswith('DetectorPlainMomentParameters'):
  s=s.replace('3/4<weighted_kappaPlain D','3/4≤weighted_kappaPlain D')
  s=s.replace('  have hΔ:0<Δ:=by linarith [D.t_pos,D.t_delta]\n','')
 if f.endswith('FinalAssemblyMomentInput'):s=s.replace(' i z Δ',' i z 0')
 if f.endswith('DetectorPlainMarkedState'):
  s=s.replace('  rw [retainedSourceEnergy_eq_energy','  change retainedSourceEnergy _ F η selected j k σ t Φ = _\n  rw [retainedSourceEnergy_eq_energy')
 if f.endswith('DetectorPlainUnmarkedField'):
  s=s.replace('  rw [actual_source_eq','  change retainedSourceEnergy (initialKeep η Q) F η ∅ j k σ t radialMajorant ≤ _\n  rw [actual_source_eq')
  s=s.replace('simpa only [hwidth,hrad,norm_zero,add_zero,one_pow,mul_one]','simpa only [hwidth,hrad,norm_zero,add_zero,one_pow,mul_one, weighted_detectorProfiles, detectorProfiles]')
 s=s.replace('(hmarked selected hselected hcap j k hjk σ hσ t ht)','(hmarked selected hselected (by simpa using hcap) j k hjk σ hσ t ht)')
 s=s.replace('≤',' ≤ ').replace('≥',' ≥ ')
 # Original modules supply all generic intermediaries and preserve actual definitions.
 s='import OAI.NumberTheory.DirichletL.'+f.replace('/','.')+'\nimport WeightedQRH.'+prev+'\n'+s
 (base/f'WeightedQRH/{target}.lean').write_text(s)
 prev=target
(base/'moment_port_order.json').write_text(json.dumps(['MomentData']+['Moment'+f.split('/')[-1] for f in files],indent=2))
print(len(names),'declarations in',len(files),'modules')
