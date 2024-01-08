'use client';
import seq from "bionode-seq";
import { useEffect, useState } from "react";
import MyButton from './MyButton';

export default function Board() {

  var [dna, setDna] = useState("");

  const [proteinSequence, setProteinSequence] = useState('');
  const [cds_seg, setCDSSeg] = useState('');
  
  const [currentTime, setCurrentTime] = useState(0);

  useEffect(() => {
    fetch('/api/time').then(res => res.json()).then(data => {
      setCurrentTime(data.time);
    });
  }, [])

  // Afegir el càlcul del percentatge GC només si s'ha introduït alguna lletra
  const gcPercentage = dna.length > 0 ? calculateGCPercentage(dna) : 0;
  console.log("Percentatge GC:", gcPercentage);

  let dna_template = seq.reverseComplement(dna);
  let rna = seq.transcribe(dna);
  let seq_prot = ""

  function handleDNA(e: Event) {
    e.preventDefault();
    let input = e.target as HTMLInputElement;
    let value = input.value.toUpperCase(); // converteix a majúscules
    let id = input.id;
    let filteredValue = value.replace(/[^ACGT]/g, '');
    console.log("filteredValue="+filteredValue)
    value = seq.reverseComplement(filteredValue);
    setDna(filteredValue);
    // Esborra la lletra d'ADN incorrecta.
    input.value = filteredValue;
  }
  
  function calculateGCPercentage(seq: string): number {
    const gcCount = (seq.match(/[GC]/g) || []).length;
    const totalBases = seq.length;
    return (gcCount / totalBases) * 100;
  }

  const handleTranslateRna = async () => {
    try {
      const response = await fetch('${process.env.NEXT_PUBLIC_API_BASE_URL}/translate', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ mrna_sequence: rna , translation_table : '1' }),
      });

      if (response.ok) {
        const data = await response.json();
        console.log(data.protein_sequence);
        setProteinSequence(data.protein_sequence);
        setCDSSeg(data.cds_seg);
      } else {
        console.error('Error en la crida POST:', response.statusText);
      }
    } catch (error) {
      console.error('Error en la crida POST:', error);
    }
  };
  
  return (
    <main className="p-6">
      <div className="bg-white shadow-md rounded px-8 pt-6 pb-8 mb-4">
        <Input id="dna-coding" pattern="/[^ACGT]/g" readOnly={false} label="DNA Coding Strand" seq={dna} handler={handleDNA} />
        <Input id="dna-template" pattern="/[^ACGT]/g" readOnly label="DNA Template strand" seq={dna_template} handler={handleDNA} />
        <Input id="rna" label="RNAm" pattern="/[^ACGU]/g" readOnly={true} seq={rna} handler={handleDNA} />
        <MyButton onClick={handleTranslateRna} />
        <p></p>
        {proteinSequence && <h2 className="block text-gray-700 text-xl pt-5 font-bold">Seqüència de proteïnes: {proteinSequence}</h2>}
        {cds_seg && <h3 className="block text-gray-700 text-xl pt-5 font-bold">Segments: {cds_seg}</h3>}
        {/* Mostra el percentatge GC només si s'ha introduït alguna lletra */}
        {dna.length > 0 && <p>Percentatge GC: {gcPercentage.toFixed(2)}%</p>}
      </div>
    </main>
  );
}


interface InputProps {
  id: string
  label: string
  seq: string
  handler: any
  readOnly: boolean
  pattern: string
}

function Input(props: InputProps) {
  return <div className="mb-4">
    <label className="block text-gray-700 text-sm font-bold mb-2">{props.label}</label>
    <input id={props.id} type="text" readOnly={props.readOnly} pattern={props.pattern} defaultValue={props.seq} 
      onKeyUp={(e) => props.handler(e)} 
      className="shadow appearance-none border rounded w-full py-2 px-3 text-gray-700 leading-tight focus:outline-none focus:shadow-outline" />
  </div>
}
