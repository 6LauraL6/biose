'use client';
import seq from "bionode-seq";
import { useEffect, useState } from "react";


export default function Board() {

  const [dna, setDna] = useState("");

  const [currentTime, setCurrentTime] = useState(0);

  useEffect(() => {
    fetch('/api/time').then(res => res.json()).then(data => {
      setCurrentTime(data.time);
    });
  }, [])



  let dna_template = seq.reverseComplement(dna);
  let rna = seq.transcribe(dna);

  function handle(e: Event) {
    e.preventDefault();
    let input = e.target as HTMLInputElement;
    let value = input.value;
    let id = input.id;
    if (id == "dna-template") {
      value = seq.reverseComplement(value);
    } else if (id == "rna") {
      // TODO
    }
    setDna(value);
  }



  return (
    <main className="p-6">
      <div className="bg-white shadow-md rounded px-8 pt-6 pb-8 mb-4">
        <Input id="dna-coding" label="DNA Coding Strand" seq={dna} handler={handle} />
        <Input id="dna-template" label="DNA Template strand" seq={dna_template} handler={handle} />
        <Input id="rna" label="RNAm" seq={rna} handler={handle} />
      </div>
      <p className="text-center mt-10">The current time is {currentTime}.</p>
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
    <input id={props.id} type="text" defaultValue={props.seq} onKeyUp={(e) => props.handler(e)} className="shadow appearance-none border rounded w-full py-2 px-3 text-gray-700 leading-tight focus:outline-none focus:shadow-outline" />
  </div>
}





