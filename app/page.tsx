'use client';
import seq from "bionode-seq";
import { useEffect, useState } from "react";


export default function Board() {

  var [dna, setDna] = useState("");

  const [currentTime, setCurrentTime,] = useState(0);

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

  function handle(e: Event) {
    e.preventDefault();
    let input = e.target as HTMLInputElement;
    let value = input.value.toUpperCase(); // converteix a majúscules
    let id = input.id;
    const filteredValue = value.replace(/[^ACGT]/g, '');
    console.log("filteredValue="+filteredValue)
    if (id == "dna-template") {
      value = seq.reverseComplement(filteredValue);
    } else if (id == "rna") {
      // TODO
    }
    setDna(filteredValue);
    // Esborra la lletra d'ADN incorrecta.
    input.value = filteredValue;
  }
  
  function calculateGCPercentage(seq: string): number {
    const gcCount = (seq.match(/[GC]/g) || []).length;
    const totalBases = seq.length;
    return (gcCount / totalBases) * 100;
  }
  
  return (
    <main className="p-6">
      <div className="bg-white shadow-md rounded px-8 pt-6 pb-8 mb-4">
        <Input id="dna-coding" label="DNA Coding Strand" seq={dna} handler={handle} />
        <Input id="dna-template" label="DNA Template strand" seq={dna_template} handler={handle} />
        <Input id="rna" label="RNAm" seq={rna} handler={handle} />
      </div>
      {/* Mostra el percentatge GC només si s'ha introduït alguna lletra */}
      {dna.length > 0 && <p>Percentatge GC: {gcPercentage.toFixed(2)}%</p>}
      <p className="text-center mt-10">The current time is {currentTime}.</p>
      <a href="./arnm-translation">Translate ARN to Proteins</a>
    </main>
  );
}


interface InputProps {
  id: string
  label: string
  seq: string
  handler: any
}

function Input(props: InputProps) {
  return <div className="mb-4">
    <label className="block text-gray-700 text-sm font-bold mb-2">{props.label}</label>
    <input id={props.id} type="text" pattern="/[^ACGT]/g" defaultValue={props.seq} onKeyUp={(e) => props.handler(e)} className="shadow appearance-none border rounded w-full py-2 px-3 text-gray-700 leading-tight focus:outline-none focus:shadow-outline" />
  </div>
}