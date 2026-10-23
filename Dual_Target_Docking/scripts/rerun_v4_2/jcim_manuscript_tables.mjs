/** Update manuscript workbooks from accepted current tables, without statistics. */
import fs from 'node:fs/promises';
import path from 'node:path';
import {createRequire} from 'node:module';
const require=createRequire(import.meta.url);
const runtime=process.env.CODEX_PRIMARY_RUNTIME_NODE_MODULES;
const entry=require.resolve('@oai/artifact-tool',runtime ? {paths:[runtime]} : undefined);
const {FileBlob,SpreadsheetFile,Workbook}=await import(entry);
const payload=JSON.parse(await fs.readFile(process.argv[2],'utf8'));
for(const spec of payload.workbooks){
 const filename=path.join(payload.root,spec.file);
 let wb;
 try{await fs.access(filename);wb=await SpreadsheetFile.importXlsx(await FileBlob.load(filename));}
 catch(e){if(e.code!=='ENOENT')throw e;wb=Workbook.create();}
 for(const item of spec.sheets){
  let sheet;
  try{sheet=wb.worksheets.getItem(item.name);}catch{sheet=wb.worksheets.add(item.name);}
  const old=sheet.getUsedRange();if(old)old.clear({applyTo:'contents'});
  const data=item.matrix;
  sheet.getRangeByIndexes(0,0,data.length,data[0].length).values=data;
  // Preserve the source export's flat table structure and fonts. Extend header
  // formatting and navigation for the new strict-subset sheet only.
  const header=sheet.getRangeByIndexes(0,0,1,data[0].length);
  if(item.name==='Table1_Label_Counts')sheet.getRange('B1:B9').format.columnWidth=21;
  if(item.name==='S3_PRIMARY_Metrics'){sheet.getRange('C1:C121').format.columnWidth=25;sheet.getRange('D1:D121').format.columnWidth=23;}
  header.format.font={name:'Arial',size:9,bold:true};
  if(item.name==='S18_Q110_Strict_Subset'){
   sheet.getRangeByIndexes(0,0,data.length,data[0].length).format.font={name:'Arial',size:9};
   header.format.font={name:'Arial',size:9,bold:true};
   sheet.getRangeByIndexes(0,0,data.length,data[0].length).format.columnWidth=19;
   sheet.freezePanes.freezeRows(1);
  }
 }
 wb.recalculate();
 console.log(spec.file,(await wb.inspect({kind:'sheet',include:'id,name',maxChars:1000})).ndjson);
 const xlsx=await SpreadsheetFile.exportXlsx(wb);await xlsx.save(filename);
 const sheetName=spec.file==='Main_Tables.xlsx'?'Table1_Label_Counts':'S3_PRIMARY_Metrics';
 const image=await wb.render({sheetName,range:spec.file==='Main_Tables.xlsx'?'A1:J9':'A1:K9',scale:1,format:'png'});
 await fs.writeFile(path.join(payload.root,'source_views',spec.file.replace('.xlsx','_preview.png')),new Uint8Array(await image.arrayBuffer()));
}
