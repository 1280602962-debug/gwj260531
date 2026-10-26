/** Publication tables only. Statistical results are copied, never recalculated. */
import fs from 'node:fs/promises';
import path from 'node:path';
import {createRequire} from 'node:module';
const require=createRequire(import.meta.url);
const runtime=process.env.CODEX_PRIMARY_RUNTIME_NODE_MODULES;
const {Workbook,SpreadsheetFile}=await import(require.resolve('@oai/artifact-tool',runtime?{paths:[runtime]}:undefined));
const payload=JSON.parse(await fs.readFile(process.argv[2],'utf8'));
await fs.mkdir(payload.preview_dir,{recursive:true});
const letters=n=>{let s='';while(n){n--;s=String.fromCharCode(65+n%26)+s;n=Math.floor(n/26);}return s;};
const textLines=text=>{const words=text.split(/\s+/);const lines=[];let line='';for(const w of words){if(line.length+w.length>145){lines.push(line);line='';}line+=(line?' ':'')+w;}if(line)lines.push(line);return lines;};
for(const kind of ['main','si']){
 const wb=Workbook.create();
 const tables=kind==='main'?payload.tables.slice(0,2):payload.tables.slice(2);
 for(const t of tables){
  const sheet=wb.worksheets.add('Table_'+t.number);sheet.showGridLines=false;
  sheet.tabColor=kind==='main'?'#284B63':'#627D98';
  const width=Math.max(...t.panels.map(p=>p.headers.length));
  sheet.getRange('A2').values=[[`Table ${t.number}. ${t.title}`]];
  sheet.getRange('A2').format.font={name:'Arial',size:14,bold:true,color:'#182B3A'};
  let row=3;
  for(const p of t.panels){
   if(p.label){sheet.getRangeByIndexes(row,0,1,1).values=[[`${p.label}. ${p.title}`]];sheet.getRangeByIndexes(row,0,1,p.headers.length).format.font={name:'Arial',size:10,bold:true};row++;}
   const matrix=[p.headers,...p.rows].map(r=>r.map(v=>typeof v==='string'&&v.startsWith('=')?"'"+v:v));
   const body=sheet.getRangeByIndexes(row,0,matrix.length,p.headers.length);body.values=matrix;
   body.format.font={name:'Arial',size:10,color:'#172B3A'};body.format.rowHeight=24;body.format.verticalAlignment='center';
   const header=sheet.getRangeByIndexes(row,0,1,p.headers.length);
   header.format.fill='#284B63';header.format.font={name:'Arial',size:10,bold:true,color:'#FFFFFF'};
   header.format.wrapText=true;header.format.rowHeight=56;header.format.horizontalAlignment='center';
   for(let c=0;c<p.headers.length;c++){
    const label=p.headers[c];
    let cw=/Target pair/.test(label)?21:/Target|PDB|Scheme|Contrast|Replaced side|Role|Metric|Cognate CCD/.test(label)?13:/Auth chains|Cognate instance/.test(label)?16:/Concordance|comparison/.test(label)?26:/condition|record|status|quota|components/i.test(label)?49:14;
    if(t.number==='S1'&&p.label==='A')cw=[21,20,56,55,15,44][c];
    sheet.getRangeByIndexes(row,c,matrix.length,1).format.columnWidth=Math.max(cw,Number(sheet.getRangeByIndexes(row,c,1,1).format.columnWidth)||0);
    const col=sheet.getRangeByIndexes(row+1,c,p.rows.length,1);
    const values=p.rows.map(r=>r[c]);
    if(values.every(v=>v===null||typeof v==='number')){
     const fraction=/fraction|recall/i.test(label);
     const integer=values.every(v=>v===null||Number.isInteger(v));
     col.setNumberFormat(fraction?'0.0%':integer?'0':'0.0000');col.format.horizontalAlignment='right';
    }else{col.format.horizontalAlignment='left';col.format.wrapText=true;}
   }
   if(t.number==='S1'&&p.label==='A'){
    sheet.getRangeByIndexes(row+1,0,p.rows.length,p.headers.length).format.rowHeight=64;
    sheet.getRangeByIndexes(row+1,1,p.rows.length,1).setNumberFormat('0.0');
   }
   if(t.number==='S2')sheet.getRangeByIndexes(row+1,0,p.rows.length,p.headers.length).format.rowHeight=55;
   if(t.number==='S10'&&p.label==='B')sheet.getRangeByIndexes(row+1,0,p.rows.length,p.headers.length).format.rowHeight=36;
   row+=p.rows.length+3;
  }
  for(const line of textLines(t.note)){sheet.getRangeByIndexes(row,0,1,1).values=[[line]];sheet.getRangeByIndexes(row,0,1,1).format.font={name:'Arial',size:10,color:'#44546A'};sheet.getRangeByIndexes(row,0,1,1).format.rowHeight=21;row++;}
  sheet.freezePanes.freezeRows(t.panels[0].label?5:4);
  sheet.freezePanes.freezeColumns(2);
 }
 wb.recalculate();
 console.log(kind,(await wb.inspect({kind:'sheet',include:'id,name',maxChars:1800})).ndjson);
 const errors=await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:30},summary:'publication table error scan',maxChars:1200});
 console.log(errors.ndjson);
 for(const t of tables){
  const cols=Math.max(...t.panels.map(p=>p.headers.length));
  const img=await wb.render({sheetName:'Table_'+t.number,range:`A1:${letters(cols)}${t.number==='S1'?14:12}`,scale:1,format:'png'});
  await fs.writeFile(path.join(payload.preview_dir,'Table_'+t.number+'.png'),new Uint8Array(await img.arrayBuffer()));
  let start=4;
  for(const [i,p] of t.panels.entries()){
   if(i>0){
    const last=letters(p.headers.length);
    const img=await wb.render({sheetName:'Table_'+t.number,range:`A${start}:${last}${start+Math.min(p.rows.length,8)+1}`,scale:1,format:'png'});
    await fs.writeFile(path.join(payload.preview_dir,'Table_'+t.number+'_'+p.label+'.png'),new Uint8Array(await img.arrayBuffer()));
   }
   start+=p.rows.length+3+(p.label?1:0);
  }
 }
 const filename=kind==='main'?'01_Main_Tables/Main_Tables.xlsx':'02_SI_Tables/Supplementary_Tables.xlsx';
 await (await SpreadsheetFile.exportXlsx(wb)).save(path.join(payload.output,filename));
 await fs.rm(path.join(payload.output,filename+'.inspect.ndjson'),{force:true});
}
