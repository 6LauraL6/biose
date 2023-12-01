'use client';
import { useEffect, useState } from "react";

export default function Board() {

  const [currentTime, setCurrentTime] = useState(0);

  useEffect(() => {
    fetch('/api/time').then(res => res.json()).then(data => {
      setCurrentTime(data.time);
    });
  }, [])


  return (
    <div className="container mt-5">
      <p className="text-center mt-5">The current time is {currentTime}.</p>
    </div>
  );
}





